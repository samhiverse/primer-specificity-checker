from flask import Flask, render_template, request, jsonify
import urllib.request
import urllib.parse
import urllib.error
import re
import time
import threading
import uuid
import html

app = Flask(__name__)

# Temporary job storage. Suitable for a small demonstration with one server worker.
blast_jobs = {}
job_lock = threading.Lock()

BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
USER_AGENT = "PrimerSpecificityChecker/1.0"
MAX_STATUS_CHECKS = 10
MIN_STATUS_INTERVAL = 60
MAX_STATUS_FAILURES = 3


def update_job(job_id, **values):
    with job_lock:
        if job_id in blast_jobs:
            blast_jobs[job_id].update(values)


def clean_dna(sequence):
    lines = [line for line in sequence.splitlines()
             if not line.strip().startswith(">")]
    return "".join("".join(lines).split()).upper()


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ATGC", "TACG"))[::-1]


def check_primers(forward, reverse, target):
    forward = clean_dna(forward)
    reverse = clean_dna(reverse)
    target = clean_dna(target)

    if not forward or not reverse or not target:
        return {"error": "Please enter all three DNA sequences."}
    if not set(forward).issubset(set("ATGC")):
        return {"error": "Forward primer must contain only A, T, G and C."}
    if not set(reverse).issubset(set("ATGC")):
        return {"error": "Reverse primer must contain only A, T, G and C."}
    if not set(target).issubset(set("ATGCN")):
        return {"error": "Target DNA may contain only A, T, G, C and N."}

    reverse_site = reverse_complement(reverse)
    f_position = target.find(forward)
    r_position = target.find(reverse_site)
    f_match = f_position != -1
    r_match = r_position != -1

    if not f_match or not r_match:
        return {
            "message": "One or both primers were not found in the target.",
            "forward_match": f_match,
            "reverse_match": r_match
        }

    if f_position < r_position:
        product_size = r_position + len(reverse_site) - f_position
        return {
            "message": "Potential target amplicon found!",
            "forward_match": True,
            "reverse_match": True,
            "product_size": product_size
        }

    return {
        "message": "Both primers matched, but their positions are not in the expected order.",
        "forward_match": True,
        "reverse_match": True
    }


def extract_hits(blast_data, primer_name):
    sections = re.split(r"(?m)^Query=\s*", blast_data)
    selected_section = None
    for section in sections:
        if section.strip().startswith(primer_name):
            selected_section = section
            break
    if selected_section is None:
        return []

    hits = []
    pattern = r"(?ms)^>([^\n]+)\n(.*?)(?=^>|\Z)"
    for match in re.finditer(pattern, selected_section):
        header = match.group(1).strip()
        alignment = match.group(2)
        accession = header.split()[0]
        strand_match = re.search(r"Strand=(Plus|Minus)/(Plus|Minus)", alignment)
        position_match = re.search(r"(?m)^Sbjct\s+(\d+)\s+[A-Z-]+\s+(\d+)", alignment)
        identity_match = re.search(r"Identities\s*=\s*(\d+)/(\d+)", alignment)
        if not strand_match or not position_match:
            continue
        query_strand, subject_strand = strand_match.groups()
        start = int(position_match.group(1))
        end = int(position_match.group(2))
        identical = int(identity_match.group(1)) if identity_match else 0
        length = int(identity_match.group(2)) if identity_match else 0
        hits.append({
            "accession": accession,
            "description": header,
            "query_strand": query_strand,
            "subject_strand": subject_strand,
            "start": min(start, end),
            "end": max(start, end),
            "identity": identical,
            "length": length
        })
    return hits


def analyse_blast(blast_data):
    forward_hits = extract_hits(blast_data, "Forward_Primer")
    reverse_hits = extract_hits(blast_data, "Reverse_Primer")
    forward_accessions = {hit["accession"] for hit in forward_hits}
    reverse_accessions = {hit["accession"] for hit in reverse_hits}
    common = forward_accessions & reverse_accessions
    products = []

    for accession in common:
        f_hits = [h for h in forward_hits if h["accession"] == accession]
        r_hits = [h for h in reverse_hits if h["accession"] == accession]
        for f in f_hits:
            for r in r_hits:
                forward_orientation = (
                    f["subject_strand"] == "Plus"
                    and r["subject_strand"] == "Minus"
                    and f["end"] < r["start"]
                )
                reverse_orientation = (
                    f["subject_strand"] == "Minus"
                    and r["subject_strand"] == "Plus"
                    and r["end"] < f["start"]
                )
                if not (forward_orientation or reverse_orientation):
                    continue
                product_size = max(f["end"], r["end"]) - min(f["start"], r["start"]) + 1
                products.append({
                    "accession": accession,
                    "description": f["description"],
                    "forward_position": f"{f['start']}-{f['end']}",
                    "reverse_position": f"{r['start']}-{r['end']}",
                    "product_size": product_size,
                    "forward_identity": f"{f['identity']}/{f['length']}",
                    "reverse_identity": f"{r['identity']}/{r['length']}"
                })

    return {
        "forward_count": len(forward_hits),
        "reverse_count": len(reverse_hits),
        "common_count": len(common),
        "products": products
    }


def request_blast(parameters, timeout=60, post=False):
    encoded = urllib.parse.urlencode(parameters)
    url = BLAST_URL if post else BLAST_URL + "?" + encoded
    payload = encoded.encode("utf-8") if post else None
    req = urllib.request.Request(
        url, data=payload,
        headers={"User-Agent": USER_AGENT, "Accept": "text/plain, text/html"}
    )
    with urllib.request.urlopen(req, timeout=timeout) as response:
        result = response.read().decode("utf-8", errors="replace")
        print("NCBI HTTP status:", response.status, flush=True)
    return result


def summarize_response(data):
    """Only log a short sanitized excerpt; never log full sequence queries."""
    without_tags = re.sub(r"<[^>]*>", " ", data)
    without_tags = html.unescape(without_tags)
    without_tags = re.sub(r"\s+", " ", without_tags).strip()
    return without_tags[:350] or "(empty response)"


def run_blast(job_id, forward, reverse):
    try:
        query = f">Forward_Primer\n{forward}\n>Reverse_Primer\n{reverse}"
        parameters = {
            "CMD": "Put",
            "PROGRAM": "blastn",
            "DATABASE": "nt",
            "QUERY": query,
            "SHORT_QUERY_ADJUST": "true",
            "WORD_SIZE": "7",
            "FILTER": "F",
            "EXPECT": "1000",
            "HITLIST_SIZE": "20",
            "FORMAT_TYPE": "Text",
            "TOOL": "PrimerSpecificityChecker"
        }
        submission = request_blast(parameters, post=True)
        rid_match = re.search(r"\bRID\s*=\s*([A-Za-z0-9_-]+)", submission)
        rtoe_match = re.search(r"\bRTOE\s*=\s*(\d+)", submission)
        if not rid_match:
            print("NCBI submission response:", summarize_response(submission), flush=True)
            raise RuntimeError("NCBI did not provide a BLAST request ID. Please try again later.")

        rid = rid_match.group(1)
        print("NCBI RID received:", rid, flush=True)
        update_job(job_id, rid=rid, status="waiting", message="NCBI is processing your search.")
        estimated_time = int(rtoe_match.group(1)) if rtoe_match else 60
        consecutive_failures = 0
        status = "WAITING"

        for attempt in range(MAX_STATUS_CHECKS):
            time.sleep(max(MIN_STATUS_INTERVAL, min(estimated_time, 120)))
            status_parameters = {
                "CMD": "Get", "RID": rid, "FORMAT_OBJECT": "SearchInfo"
            }
            try:
                status_data = request_blast(status_parameters)
                status_match = re.search(r"\bStatus\s*=\s*(WAITING|READY|FAILED|UNKNOWN)\b", status_data, re.I)
                if not status_match:
                    consecutive_failures += 1
                    print(
                        f"Unexpected NCBI status response ({consecutive_failures}/{MAX_STATUS_FAILURES}):",
                        summarize_response(status_data), flush=True
                    )
                    update_job(job_id, message="NCBI returned an unexpected status; retrying shortly.")
                    if consecutive_failures >= MAX_STATUS_FAILURES:
                        raise RuntimeError("NCBI repeatedly returned an unreadable status response. Please try again later.")
                    estimated_time = 60
                    continue
            except (urllib.error.URLError, TimeoutError, OSError) as error:
                consecutive_failures += 1
                print(f"NCBI status network error ({consecutive_failures}/{MAX_STATUS_FAILURES}): {error}", flush=True)
                update_job(job_id, message="Temporary NCBI connection problem; retrying.")
                if consecutive_failures >= MAX_STATUS_FAILURES:
                    raise RuntimeError("Unable to contact NCBI reliably. Please try again later.") from error
                estimated_time = 60
                continue

            consecutive_failures = 0
            status = status_match.group(1).upper()
            print(f"NCBI BLAST status (check {attempt + 1}): {status}", flush=True)
            update_job(job_id, message=f"NCBI search status: {status}")
            if status == "READY":
                break
            if status in ("FAILED", "UNKNOWN"):
                raise RuntimeError(f"NCBI BLAST search status: {status}")
            estimated_time = 60

        if status != "READY":
            raise RuntimeError("BLAST has not finished within the checking window. Please try again later.")

        result_parameters = {
            "CMD": "Get", "RID": rid, "FORMAT_TYPE": "Text",
            "ALIGNMENTS": "20", "DESCRIPTIONS": "20"
        }
        blast_data = request_blast(result_parameters, timeout=120)
        if not blast_data.strip() or re.search(r"\bStatus\s*=\s*WAITING\b", blast_data, re.I):
            raise RuntimeError("NCBI has not returned the completed BLAST report yet.")

        # Preserve the original report generation behavior.
        with open("ncbi_blast_results.txt", "w") as file:
            file.write(blast_data)

        analysis = analyse_blast(blast_data)
        report = [
            "PRIMER CROSS-REACTIVITY ANALYSIS",
            f"Forward primer hits: {analysis['forward_count']}",
            f"Reverse primer hits: {analysis['reverse_count']}",
            f"Common accessions: {analysis['common_count']}",
            "", "POTENTIAL PCR PRODUCT CANDIDATES"
        ]
        for number, product in enumerate(analysis["products"], start=1):
            report.extend([
                "", f"Candidate {number}",
                f"Accession: {product['accession']}",
                f"Sequence: {product['description']}",
                f"Forward position: {product['forward_position']}",
                f"Reverse position: {product['reverse_position']}",
                f"Forward identity: {product['forward_identity']}",
                f"Reverse identity: {product['reverse_identity']}",
                f"Potential product size: {product['product_size']} bp"
            ])
        if not analysis["products"]:
            report.append("No suitable paired primer arrangement was found among the parsed BLAST hits.")
        report.extend([
            "", "NOTE: Preliminary computational analysis only.",
            "This does not prove or exclude experimental amplification."
        ])
        with open("cross_reactivity_report.txt", "w") as file:
            file.write("\n".join(report))

        update_job(job_id, status="complete", message="BLAST analysis completed.", analysis=analysis)
        print("BLAST analysis completed successfully.", flush=True)

    except Exception as error:
        print("BLAST ERROR:", str(error), flush=True)
        update_job(job_id, status="error", message=str(error))


@app.route("/", methods=["GET", "POST"])
def home():
    result = None
    if request.method == "POST":
        forward = request.form.get("forward", "")
        reverse = request.form.get("reverse", "")
        target = request.form.get("target", "")
        result = check_primers(forward, reverse, target)
    return render_template("index.html", result=result)


@app.route("/start-blast", methods=["POST"])
def start_blast():
    forward = clean_dna(request.form.get("forward", ""))
    reverse = clean_dna(request.form.get("reverse", ""))
    if (not forward or not reverse
            or not set(forward).issubset(set("ATGC"))
            or not set(reverse).issubset(set("ATGC"))):
        return jsonify({
            "error": "Please enter valid forward and reverse primers using A, T, G and C."
        }), 400

    job_id = str(uuid.uuid4())
    with job_lock:
        blast_jobs[job_id] = {
            "status": "starting", "message": "Submitting your primers to NCBI."
        }
    thread = threading.Thread(target=run_blast, args=(job_id, forward, reverse), daemon=True)
    thread.start()
    return jsonify({"job_id": job_id})


@app.route("/blast-status/<job_id>")
def blast_status(job_id):
    with job_lock:
        job = blast_jobs.get(job_id)
        job = dict(job) if job is not None else None
    if not job:
        return jsonify({"error": "Search not found. Please start a new search."}), 404
    return jsonify(job)


if __name__ == "__main__":
    app.run(debug=True)


from flask import Flask, render_template, request, jsonify
import urllib.request
import urllib.parse
import re
import time
import threading
import uuid

app = Flask(__name__)

# Stores BLAST jobs while the application is running
blast_jobs = {}


# -------------------------------
# 1. CLEAN DNA SEQUENCES
# -------------------------------

def clean_dna(sequence):
    lines = [
        line for line in sequence.splitlines()
        if not line.strip().startswith(">")
    ]
    return "".join("".join(lines).split()).upper()


def reverse_complement(sequence):
    return sequence.translate(
        str.maketrans("ATGC", "TACG")
    )[::-1]


# -------------------------------
# 2. BASIC PRIMER CHECKER
# -------------------------------

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
        product_size = (
            r_position + len(reverse_site) - f_position
        )

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


# -------------------------------
# 3. EXTRACT BLAST HITS
# -------------------------------

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

        strand_match = re.search(
            r"Strand=(Plus|Minus)/(Plus|Minus)",
            alignment
        )

        position_match = re.search(
            r"(?m)^Sbjct\s+(\d+)\s+[A-Z-]+\s+(\d+)",
            alignment
        )

        identity_match = re.search(
            r"Identities\s*=\s*(\d+)/(\d+)",
            alignment
        )

        if not strand_match or not position_match:
            continue

        query_strand, subject_strand = strand_match.groups()

        start = int(position_match.group(1))
        end = int(position_match.group(2))

        if identity_match:
            identical = int(identity_match.group(1))
            length = int(identity_match.group(2))
        else:
            identical = 0
            length = 0

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


# -------------------------------
# 4. COMPARE BOTH PRIMER SEARCHES
# -------------------------------

def analyse_blast(blast_data):

    forward_hits = extract_hits(
        blast_data, "Forward_Primer"
    )

    reverse_hits = extract_hits(
        blast_data, "Reverse_Primer"
    )

    forward_accessions = {
        hit["accession"] for hit in forward_hits
    }

    reverse_accessions = {
        hit["accession"] for hit in reverse_hits
    }

    common = forward_accessions & reverse_accessions

    products = []

    for accession in common:

        f_hits = [
            h for h in forward_hits
            if h["accession"] == accession
        ]

        r_hits = [
            h for h in reverse_hits
            if h["accession"] == accession
        ]

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

                if not (
                    forward_orientation or reverse_orientation
                ):
                    continue

                product_size = (
                    max(f["end"], r["end"])
                    - min(f["start"], r["start"])
                    + 1
                )

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


# -------------------------------
# 5. SUBMIT AND RUN NCBI BLAST
# -------------------------------

def run_blast(job_id, forward, reverse):

    try:

        blast_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

        query = (
            ">Forward_Primer\n" + forward +
            "\n>Reverse_Primer\n" + reverse
        )

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

        data = urllib.parse.urlencode(
            parameters
        ).encode("utf-8")

        req = urllib.request.Request(
            blast_url,
            data=data,
            headers={
                "User-Agent": "PrimerSpecificityChecker/1.0"
            }
        )

        with urllib.request.urlopen(req, timeout=60) as response:
            submission = response.read().decode(
                "utf-8", errors="replace"
            )

            # DIAGNOSTIC LOGGING
            # This helps us understand why NCBI
            # is not returning a BLAST request ID.
            print(
                "NCBI HTTP STATUS:",
                response.status,
                flush=True
            )

            print(
                "NCBI CONTENT TYPE:",
                response.headers.get("Content-Type"),
                flush=True
            )

        print(
            "========== NCBI BLAST RESPONSE START ==========",
            flush=True
        )

        print(
            submission[:2000],
            flush=True
        )

        print(
            "========== NCBI BLAST RESPONSE END ==========",
            flush=True
        )

        rid_match = re.search(r"RID\s*=\s*(\S+)", submission)
        rtoe_match = re.search(r"RTOE\s*=\s*(\d+)", submission)

        if not rid_match:
            print(
                "ERROR: No BLAST RID found in NCBI response.",
                flush=True
            )
            raise Exception(
                "NCBI did not return a BLAST request ID."
            )

        rid = rid_match.group(1)

        print(
            "NCBI BLAST RID received:",
            rid,
            flush=True
        )

        blast_jobs[job_id]["rid"] = rid
        blast_jobs[job_id]["status"] = "waiting"
        blast_jobs[job_id]["message"] = (
            "NCBI is processing your search."
        )

        estimated_time = (
            int(rtoe_match.group(1))
            if rtoe_match else 60
        )

        # Check BLAST status no more frequently
        # than once per minute
        status = "WAITING"

        for attempt in range(10):

            time.sleep(max(60, estimated_time))

            status_parameters = {
                "CMD": "Get",
                "RID": rid,
                "FORMAT_OBJECT": "SearchInfo"
            }

            status_url = (
                blast_url + "?" +
                urllib.parse.urlencode(status_parameters)
            )

            status_req = urllib.request.Request(
                status_url,
                headers={
                    "User-Agent": "PrimerSpecificityChecker/1.0"
                }
            )

            with urllib.request.urlopen(
                status_req, timeout=60
            ) as response:
                status_data = response.read().decode("utf-8")

            status_match = re.search(
                r"Status\s*=\s*(\w+)",
                status_data
            )

            if not status_match:
                raise Exception(
                    "Could not read NCBI search status."
                )

            status = status_match.group(1)

            print(
                "NCBI BLAST STATUS:",
                status,
                flush=True
            )

            blast_jobs[job_id]["message"] = (
                f"NCBI search status: {status}"
            )

            if status == "READY":
                break

            if status in ("FAILED", "UNKNOWN"):
                raise Exception(
                    f"NCBI BLAST search status: {status}"
                )

            estimated_time = 60

        if status != "READY":
            raise Exception(
                "BLAST is taking longer than expected. "
                "Please try again later."
            )

        # Retrieve completed BLAST report
        result_parameters = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_TYPE": "Text",
            "ALIGNMENTS": "20",
            "DESCRIPTIONS": "20"
        }

        result_url = (
            blast_url + "?" +
            urllib.parse.urlencode(result_parameters)
        )

        result_req = urllib.request.Request(
            result_url,
            headers={
                "User-Agent": "PrimerSpecificityChecker/1.0"
            }
        )

        with urllib.request.urlopen(
            result_req, timeout=120
        ) as response:
            blast_data = response.read().decode("utf-8")

        # Save the BLAST report
        with open("ncbi_blast_results.txt", "w") as file:
            file.write(blast_data)

        analysis = analyse_blast(blast_data)

        # Save a simple analysis report
        report = [
            "PRIMER CROSS-REACTIVITY ANALYSIS",
            f"Forward primer hits: {analysis['forward_count']}",
            f"Reverse primer hits: {analysis['reverse_count']}",
            f"Common accessions: {analysis['common_count']}",
            "",
            "POTENTIAL PCR PRODUCT CANDIDATES"
        ]

        for number, product in enumerate(
            analysis["products"], start=1
        ):
            report.extend([
                "",
                f"Candidate {number}",
                f"Accession: {product['accession']}",
                f"Sequence: {product['description']}",
                f"Forward position: {product['forward_position']}",
                f"Reverse position: {product['reverse_position']}",
                f"Forward identity: {product['forward_identity']}",
                f"Reverse identity: {product['reverse_identity']}",
                f"Potential product size: {product['product_size']} bp"
            ])

        if not analysis["products"]:
            report.append(
                "No suitable paired primer arrangement was found "
                "among the parsed BLAST hits."
            )

        report.extend([
            "",
            "NOTE: Preliminary computational analysis only.",
            "This does not prove or exclude experimental amplification."
        ])

        with open("cross_reactivity_report.txt", "w") as file:
            file.write("\n".join(report))

        blast_jobs[job_id]["status"] = "complete"
        blast_jobs[job_id]["message"] = (
            "BLAST analysis completed."
        )
        blast_jobs[job_id]["analysis"] = analysis

        print(
            "BLAST analysis completed successfully.",
            flush=True
        )

    except Exception as error:

        print(
            "BLAST ERROR:",
            str(error),
            flush=True
        )

        blast_jobs[job_id]["status"] = "error"
        blast_jobs[job_id]["message"] = str(error)


# -------------------------------
# 6. WEBSITE ROUTES
# -------------------------------

@app.route("/", methods=["GET", "POST"])
def home():

    result = None

    if request.method == "POST":

        forward = request.form.get("forward", "")
        reverse = request.form.get("reverse", "")
        target = request.form.get("target", "")

        result = check_primers(
            forward, reverse, target
        )

    return render_template(
        "index.html",
        result=result
    )


@app.route("/start-blast", methods=["POST"])
def start_blast():

    forward = clean_dna(
        request.form.get("forward", "")
    )

    reverse = clean_dna(
        request.form.get("reverse", "")
    )

    if (
        not forward or not reverse
        or not set(forward).issubset(set("ATGC"))
        or not set(reverse).issubset(set("ATGC"))
    ):
        return jsonify({
            "error": (
                "Please enter valid forward and reverse "
                "primers using A, T, G and C."
            )
        }), 400

    job_id = str(uuid.uuid4())

    blast_jobs[job_id] = {
        "status": "starting",
        "message": "Submitting your primers to NCBI."
    }

    thread = threading.Thread(
        target=run_blast,
        args=(job_id, forward, reverse),
        daemon=True
    )

    thread.start()

    return jsonify({"job_id": job_id})


@app.route("/blast-status/<job_id>")
def blast_status(job_id):

    job = blast_jobs.get(job_id)

    if not job:
        return jsonify({
            "error": (
                "Search not found. "
                "Please start a new search."
            )
        }), 404

    return jsonify(job)


# -------------------------------
# 7. START APPLICATION
# -------------------------------

if __name__ == "__main__":
    app.run(debug=True)

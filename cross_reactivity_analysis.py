
import re

# Read the BLAST results file

try:
    with open("ncbi_blast_results.txt", "r") as file:
        blast_data = file.read()

except FileNotFoundError:
    print("ERROR: ncbi_blast_results.txt was not found.")
    print("Keep it in the same folder as this Python file.")
    raise SystemExit


# Function to extract primer matches from BLAST report

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

    # Read each alignment block

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
            identical_bases = int(identity_match.group(1))
            alignment_length = int(identity_match.group(2))
        else:
            identical_bases = 0
            alignment_length = 0

        hits.append({
            "accession": accession,
            "description": header,
            "query_strand": query_strand,
            "subject_strand": subject_strand,
            "start": min(start, end),
            "end": max(start, end),
            "identity": identical_bases,
            "alignment_length": alignment_length
        })

    return hits


# Extract forward and reverse primer matches

forward_hits = extract_hits(
    blast_data,
    "Forward_Primer"
)

reverse_hits = extract_hits(
    blast_data,
    "Reverse_Primer"
)


print()
print("=" * 60)
print("       CROSS-REACTIVITY ANALYSIS")
print("=" * 60)

print("\nForward primer matches:", len(forward_hits))
print("Reverse primer matches:", len(reverse_hits))


# Find accessions present in both searches

forward_accessions = {
    hit["accession"] for hit in forward_hits
}

reverse_accessions = {
    hit["accession"] for hit in reverse_hits
}

common_accessions = (
    forward_accessions & reverse_accessions
)

print(
    "\nSequences containing matches for both primers:",
    len(common_accessions)
)


# Compare both primer-binding sites

potential_products = []

for accession in common_accessions:

    matching_forward = [
        hit for hit in forward_hits
        if hit["accession"] == accession
    ]

    matching_reverse = [
        hit for hit in reverse_hits
        if hit["accession"] == accession
    ]

    for forward in matching_forward:

        for reverse in matching_reverse:

            f_start = forward["start"]
            f_end = forward["end"]

            r_start = reverse["start"]
            r_end = reverse["end"]

            # Check whether primers face each other
            # on the reference sequence.

            forward_orientation = (
                forward["subject_strand"] == "Plus"
                and reverse["subject_strand"] == "Minus"
                and f_end < r_start
            )

            reverse_orientation = (
                forward["subject_strand"] == "Minus"
                and reverse["subject_strand"] == "Plus"
                and r_end < f_start
            )

            if not (
                forward_orientation
                or reverse_orientation
            ):
                continue

            # Calculate possible PCR product size

            product_size = (
                max(f_end, r_end)
                - min(f_start, r_start)
                + 1
            )

            potential_products.append({
                "accession": accession,
                "description": forward["description"],
                "forward_start": f_start,
                "forward_end": f_end,
                "reverse_start": r_start,
                "reverse_end": r_end,
                "product_size": product_size,
                "forward_identity": forward["identity"],
                "forward_length": forward["alignment_length"],
                "reverse_identity": reverse["identity"],
                "reverse_length": reverse["alignment_length"]
            })


# Display and save results

report = []

report.append("PRIMER CROSS-REACTIVITY ANALYSIS")
report.append("=" * 60)

report.append(
    f"Forward primer hits analysed: {len(forward_hits)}"
)

report.append(
    f"Reverse primer hits analysed: {len(reverse_hits)}"
)

report.append(
    f"Common accessions: {len(common_accessions)}"
)

report.append("")

if potential_products:

    report.append(
        "POTENTIAL OFF-TARGET PCR PRODUCTS"
    )

    report.append("=" * 60)

    for number, product in enumerate(
        potential_products,
        start=1
    ):

        report.append("")
        report.append(f"Candidate {number}")

        report.append(
            "Accession: " + product["accession"]
        )

        report.append(
            "Sequence: " + product["description"]
        )

        report.append(
            "Forward primer position: "
            + str(product["forward_start"])
            + "-"
            + str(product["forward_end"])
        )

        report.append(
            "Reverse primer position: "
            + str(product["reverse_start"])
            + "-"
            + str(product["reverse_end"])
        )

        report.append(
            "Forward identity: "
            + str(product["forward_identity"])
            + "/"
            + str(product["forward_length"])
        )

        report.append(
            "Reverse identity: "
            + str(product["reverse_identity"])
            + "/"
            + str(product["reverse_length"])
        )

        report.append(
            "Potential product size: "
            + str(product["product_size"])
            + " bp"
        )

        report.append(
            "Interpretation: Candidate requiring further investigation."
        )

else:

    report.append(
        "No suitable paired primer arrangement was found "
        "among the parsed BLAST hits."
    )


report.append("")
report.append("=" * 60)

report.append(
    "NOTE: This is a preliminary analysis of the "
    "reported BLAST alignments."
)

report.append(
    "The report is limited to the hits returned by BLAST "
    "and does not prove or exclude experimental amplification."
)


# Print report

for line in report:
    print(line)


# Save report

with open(
    "cross_reactivity_report.txt",
    "w"
) as file:

    file.write("\n".join(report))


print()
print("Analysis completed.")
print("Report saved as cross_reactivity_report.txt")

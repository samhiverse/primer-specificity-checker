
import urllib.parse
import urllib.request

forward = "CGACAGTCAGCCGCATCTT"
reverse = "CCCCATGGTGTCTGAGCG"

candidates = {
    "PX660077.1": (20, 38, 65, 82),
    "NG_125030.1": (340, 358, 625, 642),
    "OQ420599.1": (1, 19, 46, 63)
}

def reverse_complement(sequence):
    return sequence.translate(
        str.maketrans("ATGC", "TACG")
    )[::-1]

def get_sequence(accession):
    params = urllib.parse.urlencode({
        "db": "nuccore",
        "id": accession,
        "rettype": "fasta",
        "retmode": "text"
    })

    url = (
        "https://eutils.ncbi.nlm.nih.gov/"
        "entrez/eutils/efetch.fcgi?" + params
    )

    with urllib.request.urlopen(url, timeout=30) as response:
        fasta = response.read().decode("utf-8")

    lines = fasta.strip().splitlines()
    return "".join(lines[1:]).upper()

for accession, positions in candidates.items():
    print("\nChecking:", accession)

    try:
        sequence = get_sequence(accession)

        f_start, f_end, r_start, r_end = positions

        forward_site = sequence[f_start - 1:f_end]
        reverse_site = sequence[r_start - 1:r_end]

        forward_match = forward_site == forward
        reverse_match = reverse_site == reverse_complement(reverse)

        print("Forward primer match:", forward_match)
        print("Reverse primer match:", reverse_match)

        if forward_match and reverse_match:
            product_size = r_end - f_start + 1
            print("Potential PCR product:", product_size, "bp")
        else:
            print("Candidate not verified")

    except Exception as error:
        print("Error:", error)

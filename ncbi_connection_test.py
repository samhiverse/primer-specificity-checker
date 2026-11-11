
import urllib.request
import urllib.parse

base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

parameters = {
    "db": "nuccore",
    "id": "NM_002046.7",
    "rettype": "fasta",
    "retmode": "text",
    "tool": "primer_specificity_checker"
}

url = base_url + "?" + urllib.parse.urlencode(parameters)

print()
print("=" * 50)
print("       NCBI CONNECTION TEST")
print("=" * 50)

try:
    request = urllib.request.Request(
        url,
        headers={
            "User-Agent": "PrimerSpecificityChecker/1.0"
        }
    )

    with urllib.request.urlopen(request, timeout=30) as response:
        sequence_data = response.read().decode("utf-8")

    if sequence_data.startswith(">"):

        print("\nNCBI connection successful!")
        print("\nRetrieved sequence information:")
        print(sequence_data[:300])

        lines = sequence_data.splitlines()

        dna_sequence = "".join(
            line.strip()
            for line in lines
            if not line.startswith(">")
        )

        print("\nSequence length:", len(dna_sequence), "bp")
        print("\nSequence successfully retrieved from NCBI.")

    else:
        print("NCBI returned an unexpected response.")
        print(sequence_data[:500])

except Exception as error:
    print("\nCould not retrieve the sequence.")
    print("Error:", error)

print("\n" + "=" * 50)

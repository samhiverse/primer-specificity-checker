
import urllib.request
import urllib.parse
import re
import time

# Our two test primers

forward_primer = "CGACAGTCAGCCGCATCTT"
reverse_primer = "CCCCATGGTGTCTGAGCG"

# Combine both primers into one FASTA query

query = (
    ">Forward_Primer\n"
    + forward_primer
    + "\n"
    + ">Reverse_Primer\n"
    + reverse_primer
)

# NCBI BLAST API address

blast_url = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

# BLAST search settings

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

print()
print("=" * 60)
print("          NCBI BLAST PRIMER SEARCH")
print("=" * 60)

print("\nSubmitting primers to NCBI...")

try:

    # Submit the BLAST search

    data = urllib.parse.urlencode(
        parameters
    ).encode("utf-8")

    request = urllib.request.Request(
        blast_url,
        data=data,
        headers={
            "User-Agent": "PrimerSpecificityChecker/1.0"
        }
    )

    with urllib.request.urlopen(
        request,
        timeout=60
    ) as response:

        submission = response.read().decode("utf-8")

    # Extract BLAST request ID

    rid_match = re.search(
        r"RID\s*=\s*(\S+)",
        submission
    )

    rtoe_match = re.search(
        r"RTOE\s*=\s*(\d+)",
        submission
    )

    if not rid_match:

        print("\nCould not obtain a BLAST request ID.")
        print(submission[:1000])
        raise SystemExit

    rid = rid_match.group(1)

    estimated_time = (
        int(rtoe_match.group(1))
        if rtoe_match
        else 60
    )

    print("\nBLAST search submitted successfully.")
    print("Request ID:", rid)

    print("\nWaiting for NCBI to process the search...")

    # Check search status.
    # NCBI asks users not to poll more often than once per minute.

    status = "WAITING"

    for attempt in range(10):

        wait_time = max(60, estimated_time)

        print(
            "\nChecking search status in",
            wait_time,
            "seconds..."
        )

        time.sleep(wait_time)

        status_parameters = {
            "CMD": "Get",
            "RID": rid,
            "FORMAT_OBJECT": "SearchInfo"
        }

        status_url = (
            blast_url
            + "?"
            + urllib.parse.urlencode(status_parameters)
        )

        status_request = urllib.request.Request(
            status_url,
            headers={
                "User-Agent": "PrimerSpecificityChecker/1.0"
            }
        )

        with urllib.request.urlopen(
            status_request,
            timeout=60
        ) as response:

            status_data = response.read().decode("utf-8")

        status_match = re.search(
            r"Status\s*=\s*(\w+)",
            status_data
        )

        if not status_match:

            print("Could not read the BLAST search status.")
            print(status_data[:1000])
            raise SystemExit

        status = status_match.group(1)

        print("Current status:", status)

        if status == "READY":
            break

        if status == "FAILED":

            print("\nNCBI BLAST search failed.")
            raise SystemExit

        if status == "UNKNOWN":

            print("\nNCBI could not find this request ID.")
            raise SystemExit

        if status != "WAITING":

            print("\nUnexpected BLAST status.")
            raise SystemExit

        # Subsequent checks are spaced at least one minute apart
        estimated_time = 60

    if status != "READY":

        print("\nSearch is still processing.")
        print("Run the results retrieval again later using RID:", rid)
        raise SystemExit

    # Retrieve the completed BLAST report

    print("\nRetrieving BLAST results...")

    result_parameters = {
        "CMD": "Get",
        "RID": rid,
        "FORMAT_TYPE": "Text",
        "ALIGNMENTS": "20",
        "DESCRIPTIONS": "20"
    }

    result_url = (
        blast_url
        + "?"
        + urllib.parse.urlencode(result_parameters)
    )

    result_request = urllib.request.Request(
        result_url,
        headers={
            "User-Agent": "PrimerSpecificityChecker/1.0"
        }
    )

    with urllib.request.urlopen(
        result_request,
        timeout=120
    ) as response:

        blast_results = response.read().decode("utf-8")

    # Save results to a text file

    with open(
        "ncbi_blast_results.txt",
        "w"
    ) as file:

        file.write(blast_results)

    print("\n" + "=" * 60)
    print("BLAST SEARCH COMPLETED")
    print("=" * 60)

    print("\nResults saved in ncbi_blast_results.txt")

    print("\nPreview of results:\n")

    print(blast_results[:5000])

except Exception as error:

    print("\nAn error occurred:")
    print(error)

print("\n" + "=" * 60)
print("ANALYSIS FINISHED")
print("=" * 60)

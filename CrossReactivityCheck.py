
def clean_sequence(sequence):

    return (
        sequence
        .replace(" ", "")
        .replace("\n", "")
        .replace("\r", "")
        .upper()
    )

def is_valid_dna(sequence):

    valid_bases = set("ATGCN")

    return all(
        base in valid_bases
        for base in sequence
    )

def reverse_complement(sequence):

    complement = {
        "A": "T",
        "T": "A",
        "G": "C",
        "C": "G",
        "N": "N"
    }

    return "".join(
        complement[base]
        for base in reversed(sequence)
    )

def find_matches(target, primer):

    positions = []

    start = 0

    while True:

        position = target.find(
            primer,
            start
        )

        if position == -1:
            break

        positions.append(position + 1)

        start = position + 1

    return positions

def find_amplicons(
    target,
    forward_primer,
    reverse_primer
):

    forward_positions = find_matches(
        target,
        forward_primer
    )

    reverse_binding_sequence = reverse_complement(
        reverse_primer
    )

    reverse_positions = find_matches(
        target,
        reverse_binding_sequence
    )

    amplicons = []

    for forward_position in forward_positions:

        for reverse_position in reverse_positions:

            # Reverse primer must be downstream
            if reverse_position <= forward_position:
                continue

            forward_index = forward_position - 1
            reverse_index = reverse_position - 1

            product_size = (
                reverse_index
                + len(reverse_primer)
                - forward_index
            )

            amplicons.append({
                "forward_position": forward_position,
                "reverse_position": reverse_position,
                "product_size": product_size
            })

    return amplicons

print()
print("=" * 60)
print("          PRIMER CHECKER")
print("=" * 60)

print()

target_name = input(
    "Enter target name: "
).strip()


try:

    with open(
        "target_sequence.txt",
        "r"
    ) as file:

        target_sequence = file.read()

except FileNotFoundError:

    print()
    print(
        "ERROR: target_sequence.txt was not found."
    )

    print(
        "Make sure it is in the same folder as this Python file."
    )

    exit()


target_sequence = clean_sequence(
    target_sequence
)

print()
print(
    "✓ Target sequence loaded successfully."
)

print(
    "✓ Target length:",
    len(target_sequence),
    "bp"
)

print()
print("-" * 60)
print("Enter the primer pair")
print("-" * 60)

forward_primer = input(
    "Enter forward primer: "
).strip()

reverse_primer = input(
    "Enter reverse primer: "
).strip()


forward_primer = clean_sequence(
    forward_primer
)

reverse_primer = clean_sequence(
    reverse_primer
)

print()
print("=" * 60)
print("VALIDATION")
print("=" * 60)


if not target_name:

    print("❌ Target name is missing.")
    exit()


if not target_sequence:

    print("❌ Target sequence is empty.")
    exit()


if not forward_primer:

    print("❌ Forward primer is missing.")
    exit()


if not reverse_primer:

    print("❌ Reverse primer is missing.")
    exit()


if not is_valid_dna(target_sequence):

    print(
        "❌ Target sequence contains invalid characters."
    )

    exit()


if not is_valid_dna(forward_primer):

    print(
        "❌ Forward primer contains invalid characters."
    )

    exit()


if not is_valid_dna(reverse_primer):

    print(
        "❌ Reverse primer contains invalid characters."
    )

    exit()


print("✓ All inputs are valid.")

print()
print("=" * 60)
print("TARGET INFORMATION")
print("=" * 60)

print()
print("Target:", target_name)
print(
    "Target sequence:",
    len(target_sequence),
    "bp"
)

print()
print("=" * 60)
print("FORWARD PRIMER CHECK")
print("=" * 60)

forward_matches = find_matches(
    target_sequence,
    forward_primer
)

print()
print("Forward primer:")
print(forward_primer)

print(
    "Length:",
    len(forward_primer),
    "bp"
)


if forward_matches:

    print()
    print("✓ MATCH FOUND")

    print(
        "Binding position(s):",
        forward_matches
    )

else:

    print()
    print(
        "❌ No exact match found."
    )

print()
print("=" * 60)
print("REVERSE PRIMER CHECK")
print("=" * 60)

reverse_binding_sequence = reverse_complement(
    reverse_primer
)

reverse_matches = find_matches(
    target_sequence,
    reverse_binding_sequence
)

print()
print("Reverse primer:")
print(reverse_primer)

print(
    "Length:",
    len(reverse_primer),
    "bp"
)

print()
print("Reverse-complement binding sequence:")
print(reverse_binding_sequence)


if reverse_matches:

    print()
    print("✓ MATCH FOUND")

    print(
        "Binding position(s):",
        reverse_matches
    )

else:

    print()
    print(
        "❌ No exact match found in the expected orientation."
    )

print()
print("=" * 60)
print("PRIMER PAIR CHECK")
print("=" * 60)

amplicons = find_amplicons(
    target_sequence,
    forward_primer,
    reverse_primer
)

if amplicons:

    print()
    print(
        "✓ Potential PCR product detected."
    )

    for number, amplicon in enumerate(
        amplicons,
        start=1
    ):

        print()
        print(
            "Potential amplicon",
            number
        )

        print(
            "Forward primer position:",
            amplicon["forward_position"]
        )

        print(
            "Reverse primer position:",
            amplicon["reverse_position"]
        )

        print(
            "Potential product size:",
            amplicon["product_size"],
            "bp"
        )

else:

    print()
    print(
        "❌ No suitable primer-pair arrangement detected."
    )


print()
print("=" * 60)
print("BASIC SPECIFICITY ASSESSMENT")
print("=" * 60)


if (
    forward_matches
    and reverse_matches
    and amplicons
):

    print()
    print("🟢 TARGET COMPATIBLE")

    print()
    print(
        "Both primers match the supplied target sequence."
    )

    print(
        "They are positioned in the expected orientation"
    )

    print(
        "for a potential PCR product."
    )

    print()
    print(
        "This confirms compatibility with the intended"
    )

    print(
        "target sequence."
    )

    print()
    print(
        "Note: Non-target sequences have not been checked"
    )

    print(
        "in this Update 2 prototype."
    )


elif forward_matches or reverse_matches:

    print()
    print("🟡 PARTIAL TARGET MATCH")

    print()
    print(
        "At least one primer matches the target,"
    )

    print(
        "but a complete primer-pair arrangement"
    )

    print(
        "was not detected."
    )


else:

    print()
    print("🔴 TARGET COMPATIBILITY NOT CONFIRMED")

    print()
    print(
        "The primer pair did not produce a complete"
    )

    print(
        "match on the supplied target sequence."
    )


print()
print("=" * 60)
print("ANALYSIS COMPLETE")
print("=" * 60)
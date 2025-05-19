import hashlib
import re
import sys

from Bio.Seq import Seq

segments_data = {
    1: ("PB2", 1945, 2432),
    2: ("PB1", 1945, 2432),
    3: ("PA", 1766, 2318),
    4: ("HA", 1401, 1840),
    5: ("NP", 1232, 1617),
    6: ("NA", 1147, 1506),
    7: ("MP", 801, 1052),
    8: ("NS", 692, 908),
}

geno_pat = re.compile(r"^H\d{1,2}N\d{1,2}$")

valid_nucs = set("ATGC")

seen = set()


def format_record(hdr, seq, seq_width=80):
    segment, organism, date, genotype, acc = hdr
    year = date.split("-")[0]
    genotype_string = f"({genotype})"
    org_fields = [part.strip() for part in organism.split("/")]
    if len(org_fields) > 1:
        if not genotype or genotype_string not in org_fields[-1] or geno_pat.match(genotype) is None:
            print(hdr)
            raise ValueError
        org_fields[-1] = org_fields[-1].replace(genotype_string, "").rstrip(")")
        if year and year not in org_fields[-1]:
            print(hdr)
            raise ValueError
        if org_fields[0][-2:] not in ["(A", " A"]:
            print(hdr)
            raise ValueError
    else:
        raise ValueError
    org_fields[0] = "A"
    strain = "/".join(
        [part.replace(" ", "_") for part in org_fields]
    )
    record = f">{segment}|{strain}|{genotype}|{acc}\n"
    for i in range(0, len(seq), seq_width):
        record += seq[i:i+seq_width] + "\n"
    return record

i = open(sys.argv[1])
outs = {
    v[0]: open(sys.argv[2] + "_" + v[0] + ".fa", "w")
    for v in segments_data.values()
}

total_seqs = 0
kept_seqs = 0
print(
    "Renaming segments and filtering for total length and only ATGC in seq: ..."
)
hdr = None
seq = ""
for line in i:
    if line[0] == ">":
        total_seqs += 1
        if hdr and (
            min_len <= len(seq) <= max_len
        ) and (
            all(c in valid_nucs for c in seq)
        ):
            start = seq.find("ATG")
            if start != -1 and Seq(seq[start:]).translate().find("*") > 100:
                sha256sum = hashlib.sha256(
                    seq.encode(encoding="ascii")
                ).hexdigest()
                if sha256sum not in seen:
                    try:
                        rec = format_record(hdr, seq)
                        outs[hdr[0]].write(rec)
                        seen.add(sha256sum)
                        kept_seqs += 1
                    except ValueError:
                        pass
        hdr = [part.strip() for part in line[1:].split("|")]
        try:
            segment, min_len, max_len = segments_data[int(hdr[0])]
        except (ValueError, KeyError):
            hdr = None
            continue
        hdr[0] = segment
        seq = ""
    else:
        seq += line.strip().upper()
# write last sequence
if hdr and (
    min_len <= len(seq) <= max_len
) and (
    all(c in valid_nucs for c in seq)
):
    kept_seqs += 1
    outs[hdr[0]].write(format_record(hdr, seq))

i.close()
for o in outs.values():
    o.close()

print("Done!")
print("Processed", total_seqs, "sequences,")
print("kept", kept_seqs)

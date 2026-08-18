#!/usr/bin/env python3

import getopt
import re
import sys

import pandas as pd


def main():
    # read command line arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], "i:")
    except getopt.GetoptError:
        print("isfinder_parse.py -i <ISfinder BLAST HTML page>")
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-i":
            isfinder_name = arg

    # check if all argumants were provided

    if not isfinder_name:
        print("ISfinder BLAST HTML file is not specified")
        sys.exit(1)

    # read ISfinder output html file
    isfinder_file = open(isfinder_name, "r")
    isfinder_content = isfinder_file.read()
    isfinder_file.close()

    # collect information from the table

    # NOTE (ticket 06): `re.DOTALL` is passed positionally here, which binds to
    # `re.sub`'s `count` parameter rather than `flags` -- this looks like a
    # pre-existing bug (the DOTALL flag is silently never applied). Flagged in
    # ticket 06's Comments rather than fixed, per that ticket's scope.
    isfinder_content = re.sub(  # noqa: B034
        "</article>.+", "", isfinder_content, re.DOTALL
    )  # remove HTML code after BLAST output

    isels_final = pd.DataFrame(
        columns=[
            "contig",
            "name",
            "family",
            "group",
            "origin",
            "score",
            "e-value",
            "start",
            "stop",
            "check",
        ]
    )  # specific inforamtion about every hit
    # Seed the accumulator list with the template itself (matching what
    # `isels_final.append(isels, ...)` did implicitly on each iteration):
    # `isels` is indexed by "name" while this template has "name" as an
    # ordinary column, so concatenating the two does not align the index
    # into that column -- the "name" column ends up all-NaN and a stray,
    # all-NaN "check" column survives (isels dropped "check" before this
    # point). That's pre-existing behaviour this rewrite preserves exactly,
    # not something introduced here.
    isels_final_frames = [isels_final]

    # Check if the Query= string exists
    query_exists = False

    for contig in isfinder_content.split("<b>Query=</b>")[1:]:
        query_exists = True
        isels_ge = pd.DataFrame(
            columns=["name", "family", "group", "origin"]
        )  # general information about IS
        isels = pd.DataFrame(
            columns=[
                "contig",
                "name",
                "family",
                "group",
                "origin",
                "score",
                "e-value",
                "start",
                "stop",
                "check",
            ]
        )  # specific inforamtion about every hit

        gname_match = re.match(" (\\S+)[^\n]*\n\nLength=(\\d+)\n", contig)
        gname = gname_match.group(1)  # contig name
        gsize = int(gname_match.group(2))  # contig size

        isels_ge_rows = []
        for isel in re.findall(
            "<a[^>]+>([^<]+)</a></td><td>([^<]*)</td><td>([^<]*)</td><td><a[^>]+>([^<]*)</a></td><td><a[^>]+>[^<]+</a></td><td>[^<]+</td></tr>",
            contig,
        ):  # IS name, family, group, genome
            isel_df = pd.DataFrame([list(isel)], columns=["name", "family", "group", "origin"])
            isels_ge_rows.append(isel_df)
        # Only concat when there is at least one row: concatenating with the
        # still-empty `isels_ge` accumulator is unnecessary and triggers
        # pandas's "empty/all-NA entries" FutureWarning; skipping it when
        # there's nothing to add keeps `isels_ge` as the empty, correctly
        # columned DataFrame it already is.
        if isels_ge_rows:
            isels_ge = pd.concat(isels_ge_rows, sort=False)

        isels_ge = isels_ge.set_index("name")

        # collect coordinates of IS elements
        match = re.search("<pre>&gt;(.+)", contig, re.DOTALL)
        if match is None:
            continue

        alignments_block = match.group(1)
        isels_rows = []
        for alignments in alignments_block.split("&gt;"):
            match = re.search(r"</a> (\S+) <span", alignments)
            name = match.group(1)
            index = 1

            for align in alignments.split(" Score")[1:]:
                evalue = float(re.search("Expect = ([^\n]+)", align).group(1))
                if evalue > 1e-30:
                    continue
                bit = int(re.search(r" = (\d+) bits", align).group(1))
                match = re.findall("Query[^\n]+\n", align)
                smatch = re.search(r"Query\s+(\d+)\s+", match[0])
                start = smatch.group(1)
                smatch = re.search(r"Query\s+\d+\s+[ATGCNatgcn-]+\s+(\d+)", match[-1])
                try:
                    stop = smatch.group(1)
                except Exception:
                    print(name)
                    print(match)

                add_list = [gname, name + "_" + str(index)]
                add_list.extend(list(isels_ge.loc[name, :]))
                add_list.extend([bit, evalue, start, stop, 0])
                isel_df = pd.DataFrame(
                    [add_list],
                    columns=[
                        "contig",
                        "name",
                        "family",
                        "group",
                        "origin",
                        "score",
                        "e-value",
                        "start",
                        "stop",
                        "check",
                    ],
                )
                isels_rows.append(isel_df)
                index = index + 1

        # Same rationale as isels_ge above: skip concatenating when nothing
        # was collected so `isels` stays the empty, correctly columned
        # DataFrame it was initialized as.
        if isels_rows:
            isels = pd.concat(isels_rows, sort=False)

        isels = isels.sort_values(by=["score"], ascending=False)
        isels = isels.set_index("name")

        genome = [0] * gsize  # set genome coverage list

        for index in isels.index:
            genome_temp = genome.copy()
            is_new = 0  # new space occupied by IS element
            stop = int(isels.at[index, "stop"])
            start = int(isels.at[index, "start"])
            is_length = stop - start  # length of the IS element
            for pos in range(start - 1, stop):
                if genome_temp[pos] == 0:
                    is_new = is_new + 1

                genome_temp[pos] = 1

            if (
                is_new / is_length > 0.75
            ):  # check if IS element occupies a new region at least 75% of its length
                genome = genome_temp
                isels.at[index, "check"] = 1

        isels = isels[isels["check"] == 1]  # remove overlapping IS hits with less score
        isels = isels.drop(columns=["check"])  # remove check column
        # Only collect non-empty per-contig tables: an empty table
        # contributes no rows, so skipping it changes nothing but the
        # (irrelevant) presence of a no-op entry in the concat list.
        if not isels.empty:
            isels_final_frames.append(isels)

    # Unlike the isels_ge/isels concats above, this one cannot skip the
    # empty-accumulator case: isels_final_frames always starts with the
    # (empty) isels_final template (see the comment where it's built), which
    # is what produces the blank "name"/stray "check" columns being
    # preserved. That means this call does trigger pandas's "empty/all-NA
    # entries" FutureWarning even when isels were found -- expected here,
    # not an oversight.
    isels_final = pd.concat(isels_final_frames, sort=False)

    if query_exists:
        # write full table in file
        isels_final.to_csv("ISTable_full.txt", sep="\t")
        # write table for further processing to file
        isels_final.to_csv(
            "ISTable_processing.txt", sep="\t", columns=["contig", "start", "stop"], header=False
        )
    else:
        print(
            "\nWARNING!:"
            'No "Query=" strings with contig names found in the provided page.\n'
            "ISFinder BLAST requires short Fasta headers to include contig names.\n"
            "Please make shorter contig names or remove auxiliary information from the"
            " Fasta header."
        )


if __name__ == "__main__":
    main()

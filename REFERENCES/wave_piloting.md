# Wave piloting (Marshall Islands)

Pointers, not results. Nothing here was computed in this repository, and
nothing here was read from its source: every page was out of reach of the
build environment on 2026-11-15. Co-authored. CC0.

Read levels, as in `related_studies.md`:

```text
ABSTRACT READ   the abstract page was fetched; the body was not read
NOT FETCHED     the page could not be reached from the build environment
CARRIED         a statement supplied by a co-author, not checked against the
                source here
```

A line under "carried" is a reason to read the source. It is not a citation
for the statement. Titles and author lists that were not supplied are left
void; they are not filled in from memory.

---

## W1. Genz et al. 2009

```text
where     Oceanography 22(2):234-245, 2009. doi 10.5670/oceanog.2009.52
status    NOT FETCHED (tos.org failed to resolve; doi.org and crossref were
          already out of reach). Title, full author list and abstract: void.
carried   navigators detect reflected swell 40 km out;
          that is below buoy detection.
```

## W2. Genz 2014

```text
where     Ethos, 2014. Volume, pages and doi: void (not supplied, not fetched).
status    NOT FETCHED (the publisher's host failed to resolve). Title and
          abstract: void.
carried   egocentric sensing + allocentric map.
```

## W3. Genz 2025 voyage study

```text
where     University of Hawaii at Hilo. No document identified.
status    NOT FETCHED (hilo.hawaii.edu failed to resolve). Whether anything
          is published: void.
carried   a 2025 voyage study; ongoing.
```

## W4. Great Lakes crest-angle reading

```text
where     no document. Oral.
status    CARRIED. There is no source to fetch.
carried   reading the angle of wave crests for the wind's trend.
          Offered as a parallel in structure only, not as the same practice
          and not as evidence for W1 - W3.
```

---

Pointer to: nothing in the toolkit yet. No run here reads a wave field.
The nearest pieces are `field/hexsheet.py` (a sheet answering a field, and
`inverse_readout`, which recovers a field from the sheet's response) and
`docs/resolution_leak.md` (a reading bounded by its instrument). Neither was
built from these sources, and no link between them and wave piloting is
claimed.

Upstream and downstream ORFs
============================

Search model
------------

``uorf`` and ``dorf`` reconstruct transcript and CDS sequences before scanning.
Initiation uses ``ATG`` and termination uses the standard ``TAA``, ``TAG`` and
``TGA`` codons in the candidate's reading frame. There is no CLI option for
near-cognate starts or an alternative mitochondrial genetic code.

uORF candidates start in the 5-prime leader. The implementation distinguishes:

* ``type1``: termination within the leader, before the canonical CDS.
* ``type2``: termination after crossing the CDS start but before the canonical
  CDS end.
* ``type3``: an in-frame upstream initiation candidate extending to the canonical
  CDS stop, interpreted as an N-terminal extension candidate.

dORFs start downstream of the canonical CDS and terminate within the 3-prime
sequence. Sequence-defined ORFs are not evidence of ribosome occupancy or
protein expression. Annotation completeness affects which candidates exist.

Run extraction
--------------

.. code-block:: bash

   GenoKit uorf -d results/human.gtf.sqlite -g human.fa -s gtf -l 6 -f csv -o results/uorf.csv
   GenoKit dorf -d results/human.gtf.sqlite -g human.fa -s gtf -l 6 -f gff -o results/dorf.gff

``-l`` is measured in nucleotides, not amino acids. The default is 6; output
keeps lengths **strictly greater than** the threshold. A 6-nt candidate is
therefore excluded with ``-l 6``. Report this rule when summarizing ORF numbers.

Coordinates and figures
-----------------------

CSV reports transcript-relative ORF locations; GFF/GTF output projects the ORF
onto genomic blocks across splice junctions. Do not count each genomic block
as a separate ORF. See :doc:`inputs_outputs` for the legacy CDS interval field.

For a selected transcript, the parsers accept string-valued schematic options:

.. code-block:: bash

   GenoKit uorf -d results/human.gtf.sqlite -g human.fa -s gtf \
     -i TRANSCRIPT_ID -l 6 -f csv -o results/selected_uorf.csv \
     -m results/uorf_spliced -n results/uorf_genomic

``-m`` and ``-n`` require output prefixes; they are not standalone Boolean
switches. Use ``GenoKit view`` for an overview of every transcript in a gene.

Why a transcript has no ORF result
----------------------------------

Possible reasons include no usable CDS, an absent leader/trailer, no ATG-stop
combination, or candidates below the threshold. The detector locates the CDS
within the reconstructed transcript sequence; incomplete, inconsistent or
repeated sequence contexts should be inspected before biological interpretation.
Do not treat missing output as proof that translation cannot occur.

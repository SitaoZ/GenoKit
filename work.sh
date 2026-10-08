# requirements
pip freeze > requirements.txt
# dist 
python -m pip install build
python -m build --sdist

# Gene ID
Extract
- Gene、Promoter、Terminator、IGR 显示对应物种的 Gene ID。
- Transcript、Exon、Intron、CDS、UTR、uORF、dORF 显示对应物种的 Transcript ID。

Design
- Primer、sgRNA、shRNA、siRNA 自动切换为对应物种的 Gene ID。
- Barcode、Motif 不受影响。

Visualize
- Gene structure 自动切换为对应物种的 Gene ID。

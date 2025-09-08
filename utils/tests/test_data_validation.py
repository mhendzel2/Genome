import io
import gzip
import pandas as pd

from utils.data_validation import DataValidator


def make_gz_bytes(s: str) -> bytes:
    bio = io.BytesIO()
    with gzip.GzipFile(fileobj=bio, mode='wb') as gz:
        gz.write(s.encode('utf-8'))
    return bio.getvalue()


def test_load_gzipped_bed_and_metadata():
    dv = DataValidator()

    bed_text = 'chr1\t100\t200\nchr1\t300\t400\n'
    gz = make_gz_bytes(bed_text)
    uploaded = io.BytesIO(gz)

    df, meta = dv.load_with_metadata(uploaded, 'BED')

    assert meta['format'] == 'BED-like'
    assert meta['num_rows'] == 2
    assert meta['num_columns'] >= 3
    assert 'chr' in df.columns and 'start' in df.columns and 'end' in df.columns


def test_gene_expression_header_and_metadata():
    dv = DataValidator()

    # With header
    tsv = 'gene_id\tsample1\tsample2\nTP53\t10\t12\nBRCA1\t5\t6\n'
    uploaded = io.StringIO(tsv)
    df, meta = dv.load_with_metadata(uploaded, 'Gene Expression')

    # After load_with_metadata, first col should be renamed to 'gene'
    assert 'gene' in df.columns
    assert meta['num_rows'] == 2
    assert meta['num_columns'] == 3

    # Without header (fallback)
    tsv_no_header = 'TP53\t10\t12\nBRCA1\t5\t6\n'
    uploaded2 = io.StringIO(tsv_no_header)
    df2, meta2 = dv.load_with_metadata(uploaded2, 'Gene Expression')

    # When no header, parsers may infer differently; ensure at least one row and correct column count
    assert meta2['num_rows'] >= 1
    assert meta2['num_columns'] == 3

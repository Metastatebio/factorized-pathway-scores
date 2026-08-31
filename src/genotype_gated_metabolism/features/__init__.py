"""Cross-omics feature and signature construction."""

from .signatures import build_expression_signature, parse_gene_signature

__all__ = ["build_expression_signature", "parse_gene_signature"]

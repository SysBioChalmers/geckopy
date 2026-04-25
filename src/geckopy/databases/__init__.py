"""Database loaders for organism-specific data used by geckopy."""
from .complex_portal_loader import ComplexPortalEntry, load_complex_portal_json
from .uniprot_loader import UniprotDB, load_uniprot_tsv

__all__ = [
    "ComplexPortalEntry",
    "UniprotDB",
    "load_complex_portal_json",
    "load_uniprot_tsv",
]

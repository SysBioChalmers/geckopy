"""Database loaders for organism-specific data used by geckopy."""
from .uniprot_loader import UniprotDB, load_uniprot_tsv

__all__ = ["UniprotDB", "load_uniprot_tsv"]

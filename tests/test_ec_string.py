"""Tests for get_ec_string."""
import logging

import pytest

from geckopy.get_enzyme_data import get_ec_string


# --------------------------------------------------------------------------- #
# Trivial cases
# --------------------------------------------------------------------------- #

def test_empty_string_returns_empty():
    assert get_ec_string("") == ""


def test_whitespace_only_returns_empty(caplog):
    with caplog.at_level(logging.WARNING):
        assert get_ec_string("   \t\n") == ""
    # All-whitespace splits to no tokens; nothing was actually skipped.
    assert "skipped" not in caplog.text


def test_single_token():
    assert get_ec_string("1.1.1.1") == "EC1.1.1.1"


def test_multiple_tokens():
    assert get_ec_string("1.1.1.1 2.7.7.7") == "EC1.1.1.1 EC2.7.7.7"


def test_three_tokens_matches_header_example():
    """The MATLAB header docs give the example
    `ECX.X.X.X ECX.X.X.X ECX.X.X.X`."""
    assert (
        get_ec_string("1.1.1.1 2.7.7.7 3.4.21.1")
        == "EC1.1.1.1 EC2.7.7.7 EC3.4.21.1"
    )


# --------------------------------------------------------------------------- #
# Whitespace handling (MATLAB strsplit collapse-delimiters parity)
# --------------------------------------------------------------------------- #

def test_runs_of_whitespace_collapsed():
    assert get_ec_string("1.1.1.1   2.7.7.7") == "EC1.1.1.1 EC2.7.7.7"


def test_leading_and_trailing_whitespace_dropped():
    assert get_ec_string("  1.1.1.1 2.7.7.7  ") == "EC1.1.1.1 EC2.7.7.7"


def test_mixed_whitespace_kinds():
    assert get_ec_string("1.1.1.1\t2.7.7.7\n3.4.21.1") == (
        "EC1.1.1.1 EC2.7.7.7 EC3.4.21.1"
    )


# --------------------------------------------------------------------------- #
# Semicolon stripping (MATLAB element-wise replace parity)
# --------------------------------------------------------------------------- #

def test_trailing_semicolon_stripped():
    assert get_ec_string("1.1.1.1;") == "EC1.1.1.1"


def test_semicolon_after_each_token_stripped():
    assert get_ec_string("1.1.1.1; 2.7.7.7;") == "EC1.1.1.1 EC2.7.7.7"


def test_lone_semicolon_skipped_with_warning(caplog):
    with caplog.at_level(logging.WARNING):
        result = get_ec_string(";")
    assert result == ""
    assert "skipped" in caplog.text


# --------------------------------------------------------------------------- #
# Already-EC-prefixed input (idempotency)
# --------------------------------------------------------------------------- #

def test_already_prefixed_token_is_idempotent():
    assert get_ec_string("EC1.1.1.1") == "EC1.1.1.1"


def test_round_trip_idempotent():
    once = get_ec_string("1.1.1.1 2.7.7.7")
    twice = get_ec_string(once)
    assert once == twice == "EC1.1.1.1 EC2.7.7.7"


def test_lowercase_ec_prefix_stripped():
    assert get_ec_string("ec1.1.1.1") == "EC1.1.1.1"


def test_mixed_case_ec_prefix_stripped():
    assert get_ec_string("Ec1.1.1.1 eC2.7.7.7") == "EC1.1.1.1 EC2.7.7.7"


def test_lone_EC_prefix_skipped_with_warning(caplog):
    """A token of just `EC` strips to `""` and should be flagged."""
    with caplog.at_level(logging.WARNING):
        result = get_ec_string("EC")
    assert result == ""
    assert "skipped" in caplog.text
    assert "'EC'" in caplog.text


# --------------------------------------------------------------------------- #
# IUBMB `-` placeholders accepted in any position
# --------------------------------------------------------------------------- #

def test_dash_in_last_level_accepted():
    assert get_ec_string("1.1.1.-") == "EC1.1.1.-"


def test_dash_cascade_accepted():
    assert get_ec_string("1.-.-.-") == "EC1.-.-.-"


def test_all_dashes_accepted():
    assert get_ec_string("-.-.-.-") == "EC-.-.-.-"


def test_mixed_real_and_dashed():
    assert get_ec_string("1.1.1.1 1.1.1.- 1.-.-.-") == (
        "EC1.1.1.1 EC1.1.1.- EC1.-.-.-"
    )


# --------------------------------------------------------------------------- #
# Multi-digit levels (e.g. EC 3.4.21.1)
# --------------------------------------------------------------------------- #

def test_multi_digit_levels():
    assert get_ec_string("3.4.21.1 2.7.11.1") == "EC3.4.21.1 EC2.7.11.1"


# --------------------------------------------------------------------------- #
# Validation: invalid tokens skipped with warning
# --------------------------------------------------------------------------- #

def test_three_level_token_rejected(caplog):
    with caplog.at_level(logging.WARNING):
        result = get_ec_string("1.1.1")
    assert result == ""
    assert "skipped" in caplog.text
    assert "'1.1.1'" in caplog.text


def test_five_level_token_rejected(caplog):
    with caplog.at_level(logging.WARNING):
        result = get_ec_string("1.1.1.1.1")
    assert result == ""
    assert "skipped" in caplog.text


def test_alphabetic_junk_rejected(caplog):
    with caplog.at_level(logging.WARNING):
        result = get_ec_string("notanec")
    assert result == ""
    assert "skipped" in caplog.text
    assert "'notanec'" in caplog.text


def test_partially_alphabetic_token_rejected(caplog):
    with caplog.at_level(logging.WARNING):
        result = get_ec_string("1.1.1.x")
    assert result == ""
    assert "skipped" in caplog.text


def test_negative_number_rejected(caplog):
    """`-1` would tempt the regex; the dash is for placeholders only."""
    with caplog.at_level(logging.WARNING):
        result = get_ec_string("1.1.1.-1")
    assert result == ""
    assert "skipped" in caplog.text


def test_mixed_valid_and_invalid_keeps_valid(caplog):
    with caplog.at_level(logging.WARNING):
        result = get_ec_string("1.1.1.1 junk 2.7.7.7")
    assert result == "EC1.1.1.1 EC2.7.7.7"
    assert "skipped" in caplog.text
    assert "'junk'" in caplog.text


def test_warning_lists_all_invalid_tokens(caplog):
    with caplog.at_level(logging.WARNING):
        get_ec_string("junk1 1.1.1.1 junk2 junk3")
    assert "junk1" in caplog.text
    assert "junk2" in caplog.text
    assert "junk3" in caplog.text


def test_no_warning_when_all_tokens_valid(caplog):
    with caplog.at_level(logging.WARNING):
        get_ec_string("1.1.1.1 2.7.7.7")
    assert "skipped" not in caplog.text


# --------------------------------------------------------------------------- #
# Combined: prefix + semicolon + whitespace + dash
# --------------------------------------------------------------------------- #

def test_realistic_uniprot_style_export():
    """A token mixing all the warts seen in real DB exports."""
    raw = "  EC1.1.1.1;   EC2.7.7.-   "
    assert get_ec_string(raw) == "EC1.1.1.1 EC2.7.7.-"

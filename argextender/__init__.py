"""ARGextender: iterative read recruitment and targeted assembly around antimicrobial resistance genes."""

__version__ = "2.0.0.dev0"


class ArgExtenderError(Exception):
    """An error with a message meant for the user: bad input, a missing tool or a failed tool call."""

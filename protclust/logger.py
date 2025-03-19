import logging

logger = logging.getLogger("protein_clustering")
handler = logging.StreamHandler()
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.info)


def set_verbosity(verbose=0):
    """Set the verbosity level of the logger.

    Args:
        verbose (int): If 0, set to DEBUG level. If 1, INFO level. Else, WARNING level.
    """
    if verbose == 2:
        logger.setLevel(logging.DEBUG)
    elif verbose == 1:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.WARNING)


def get_verbosity() -> str:
    """
    Get the current verbosity level of the logger.

    """
    return logging.getLevelName(logger.level)

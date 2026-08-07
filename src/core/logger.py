##############################################################
# Logger
##############################################################

from __future__ import annotations

import logging
from pathlib import Path
from typing import Dict, List, Literal, Optional, Union, Any


class ConditionFormatter(logging.Formatter):
    """
    Formatter for the Log Message

    Attributes
    ----------
    logging.Formatter
        Logging Formatter
    """
    def format(self, record: logging.LogRecord) -> str:
        """
        Extra Formater to add a supplementary | data_name |
        
        Parameters
        ----------
        record: logging.LogRecord
            Record of the Logger

        Returns
        -------
        str
            Logger Output
        """
        base = f"{self.formatTime(record)} | {record.levelname} | {record.name}"

        if hasattr(record, "data_name"):
            base += f" | {record.data_name}"

        base += f" | {record.getMessage()}"

        return base
    

def setup_logging(logfile: Union[str,None] = None, level: int = logging.INFO) -> None:
    """
    Setup the Logging layout

    Parameters
    ----------
    logfile: str or None, default=None
        Str to the logfile
    level: int, default=logging.INFO
        Level of the logger, e.g. info, error, warning

    Returns
    -------
    None
    """

    logger = logging.getLogger()
    logger.setLevel(level)

    if logger.handlers:
        return

    formatter = ConditionFormatter()

    # Console output
    console = logging.StreamHandler()
    console.setFormatter(formatter)
    logger.addHandler(console)

    if logfile is not None:
        root_name = logging.getLogger().name
        if root_name == "Spock_0":
            logfile = logfile + "_tmp"
        logfile = Path(logfile)
        logfile.parent.mkdir(parents=True, exist_ok=True)

        file_handler = logging.FileHandler(logfile)
        file_handler.setFormatter(formatter)
        logger.addHandler(file_handler)

        cwd_logfile = Path.cwd() / "Tiberius.log"
        cwd_handler = logging.FileHandler(cwd_logfile)
        cwd_handler.setFormatter(formatter)
        logger.addHandler(cwd_handler)
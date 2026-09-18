import logging

from rich.logging import RichHandler

SUCCESS_LEVEL_NUM = 25
logging.addLevelName(SUCCESS_LEVEL_NUM, "SUCCESS")


class Logger:
    """
    A singleton console logger for ChemSampler, following the Ersilia logging pattern.
    """

    _instance = None

    def __new__(cls, *args, **kwargs):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
        return cls._instance

    def __init__(self):
        if hasattr(self, "logger"):
            return
        self.logger = logging.getLogger("chemsampler")
        self.logger.setLevel(logging.INFO)
        self.logger.handlers.clear()
        handler = RichHandler(
            rich_tracebacks=True,
            markup=False,
            log_time_format="%H:%M:%S",
            show_path=False,
        )
        handler.setFormatter(logging.Formatter("%(message)s"))
        self.logger.addHandler(handler)

    def debug(self, text):
        self.logger.debug(text)

    def info(self, text):
        self.logger.info(text)

    def warning(self, text):
        self.logger.warning(text)

    def error(self, text):
        self.logger.error(text)

    def critical(self, text):
        self.logger.critical(text)

    def success(self, text):
        """Log a message at the SUCCESS level (between INFO and WARNING)."""
        self.logger.log(SUCCESS_LEVEL_NUM, text)


logger = Logger()

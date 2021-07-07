import logging
import subprocess
from typing import List

logger = logging.getLogger(__name__)


def execute_commands(commands: List[List[str]]):
    try:
        for command in commands:
            logger.debug(f'Running [{" ".join(command)}]')
            subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True, text=True)
    except subprocess.CalledProcessError as e:
        err_msg = str(e.stderr.strip())
        raise Exception(f'Could not run [{" ".join(e.cmd)}]: error {err_msg}')

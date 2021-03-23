from pathlib import Path
import subprocess


class VariationFile:

    def __init__(self, file: Path):
        self._file = file

    def write(self, output: Path, file_type: str = 'bcf') -> Path:
        if file_type == 'bcf':
            command_bcf = ['bcftools', 'view', str(self._file), '-o', str(output), '-O', 'b']
            command_index = ['bcftools', 'index', str(output)]
            try:
                subprocess.run(command_bcf, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                           check=True, text=True)
                subprocess.run(command_index, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                           check=True, text=True)
                return output
            except subprocess.CalledProcessError as e:
                err_msg = str(e.stderr.strip())
                raise Exception(f'Could not run [{" ".join(e.cmd)}] on original file [{self._file}]: error {err_msg}')
        else:
            raise Exception(f'Invalid file_type=[{file_type}]')

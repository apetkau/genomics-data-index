from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


class SamplesProgressLogger:
    """
    Used encapsulate progress monitoring as I'm loading in data.
    """
    total_batches = None
    current_batch = 0
    total_stages = 2

    def __init__(self, stage_name: str, stage_number: int, total_samples: int):
        self._total_samples = total_samples
        self._stage_name = stage_name
        self._stage_number = stage_number

    def get_update_every_nth_sample(self, percent_to_update: float) -> int:
        if percent_to_update < 0 or percent_to_update > 100:
            raise Exception(f'percent_to_update={percent_to_update} must be between 0 and 100')

        sample_to_update = int(self._total_samples * (percent_to_update / 100))
        return max(1, sample_to_update)

    def is_batch_mode(self) -> bool:
        return self.total_batches is not None

    def update_progress(self, sample_number: int) -> None:
        if self.is_batch_mode():
            prefix = f'Sample batch {self.current_batch}/{self.total_batches}: '
        else:
            prefix = ''
        logger.info(f'{prefix}Stage {self._stage_number}/{self.total_stages} ({self._stage_name}): '
                    f'Processed {sample_number / self._total_samples * 100:0.0f}% '
                    f'({sample_number}/{self._total_samples}) samples')

    @classmethod
    def set_total_batches(cls, total_batches: int) -> None:
        cls.total_batches = total_batches

    @classmethod
    def set_current_batch(cls, current_batch: int) -> None:
        cls.current_batch = current_batch

    @classmethod
    def unset_batch_mode(cls) -> None:
        cls.total_batches = None

from __future__ import annotations

import logging

logger = logging.getLogger(__name__)


class SamplesProgressLogger:
    """
    Used encapsulate progress monitoring as I'm loading in data.
    """

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

    def update_progress(self, sample_number: int) -> None:
        logger.info(f'Stage {self._stage_number}/{self.total_stages} ({self._stage_name}): '
                    f'Processed {sample_number / self._total_samples * 100:0.0f}% '
                    f'({sample_number}/{self._total_samples}) samples')

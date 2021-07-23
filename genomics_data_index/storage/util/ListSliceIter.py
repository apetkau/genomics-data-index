import math
from typing import List, Any, Generator


class ListSliceIter:

    def __init__(self, data: List[Any], slice_size: int = 1):
        if not isinstance(data, list):
            raise Exception(f'data={data} is not of type list')
        if not isinstance(slice_size, int) or slice_size <= 0:
            raise Exception(f'slice_size={slice_size} must be a positive integer')

        self._data = data
        self._slice_size = slice_size

    def islice(self) -> Generator[List[Any], None, None]:
        num_chunks = math.ceil(len(self._data) / self._slice_size)
        for chunk in range(num_chunks):
            start = chunk * self._slice_size
            end = min(start + self._slice_size, len(self._data))
            yield self._data[start:end]

import abc

from storage.variant.io.FeaturesReader import FeaturesReader


class SampleFilesProcessor(abc.ABC):

    def __init__(self):
        pass

    @abc.abstractmethod
    def preprocess_files(self, reader: FeaturesReader) -> FeaturesReader:
        pass

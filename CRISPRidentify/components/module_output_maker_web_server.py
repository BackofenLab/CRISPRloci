import os
from os import listdir
from os.path import isfile, join

from components_output_maker import SummaryMakerCSV
from components_output_maker import GFFOutputMaker
from components_output_maker_web_server import SimpleOutputMakerWebServer
from components_output_maker_web_server import SummaryOutputMakerWebServer
from components_output_maker_web_server import AdditionalOutputWebServer
from components_output_maker_web_server import CRISPRArraysInSingleFiles


class OutputMaker:
    def __init__(self, file_path, parameters, result_path, pickle_result_path,
                 categories, non_array_data, list_features, header):
        self.file_path = file_path
        self.parameters = parameters
        self.result_path = result_path
        self.pickle_result_path = pickle_result_path
        self.categories = categories
        self.non_array_data = non_array_data
        self.list_features = list_features
        self.header = header

        self._make_output()
        self._remove_empty_files()

    def _make_output(self):
        som = SimpleOutputMakerWebServer(categories=self.categories,
                                         result_path=self.result_path,
                                         non_array_data=self.non_array_data,
                                         list_features=self.list_features)

        suom = SummaryOutputMakerWebServer(result_path=self.result_path,
                                           categories=self.categories,
                                           non_array_data=self.non_array_data,
                                           header=self.header,
                                           list_feature_names=self.list_features)

        sm_csv = SummaryMakerCSV(result_path=self.result_path,
                                 categories=self.categories,
                                 non_array_data=self.non_array_data)

        aows = AdditionalOutputWebServer(categories=self.categories,
                                         non_array_data=self.non_array_data,
                                         result_path=self.result_path)

        gffom = GFFOutputMaker(result_path=self.result_path,
                               categories=self.categories,
                               non_array_data=self.non_array_data,
                               header=self.header,
                               list_feature_names=self.list_features)

        casf = CRISPRArraysInSingleFiles(result_path=self.result_path,
                               categories=self.categories,
                               non_array_data=self.non_array_data,
                               header=self.header,
                               list_feature_names=self.list_features)



    def _remove_empty_files(self):
        files_in_result_dir = [join(self.result_path, f) for f in listdir(self.result_path)
                               if isfile(join(self.result_path, f))]

        for file_path in files_in_result_dir:
            size = os.path.getsize(file_path)
            if size == 0:
                try:
                    os.remove(file_path)
                except Exception:
                    pass

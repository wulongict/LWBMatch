#!/usr/bin/python
# Python 2.7 Required
# Build the network of depencencies among each step in the entire workflow.
# A object named as wfc, keeps a binary variable for each step. True for run, False for skip.
#
# -----------------------
# Author: L. WU
# Date: Aug. 18, 2016
# ------------------------




class WorkflowController:
    def __init__(self):
        self.run_indexmzXML = False

        self.run_dia_umpire = False
        self.run_ms_convert = False
        self.run_db_generate = False
        self.run_comet = False
        self.run_xinteract = False
        self.run_xinteractQ1 = True
        self.run_xinteractQ2 = False
        self.run_xinteractQ3 = False
        self.run_extract_ground_truth = False

        self.run_feature_finder = False
        self.run_lwbmatch = False
        self.run_plot_figure = False

        self.run_recall_precision = False
        self.run_generate_random_feature = False

        self.run_dtw = False

    def reset_workflow(self):
        if self.run_indexmzXML:
            self.run_feature_finder = True
            self.run_dia_umpire = True

        if self.run_generate_random_feature:
            self.run_lwbmatch = True

        run_dia_umpire_flow = [self.run_dia_umpire, self.run_ms_convert, self.run_db_generate, self.run_comet,
                               self.run_xinteract, self.run_extract_ground_truth]
        run_lwbmatch_flow = [self.run_feature_finder, self.run_lwbmatch]

        for i in range(1, len(run_dia_umpire_flow)):
            if run_dia_umpire_flow[i - 1]:
                run_dia_umpire_flow[i] = True

        for i in range(1, len(run_lwbmatch_flow)):
            if run_lwbmatch_flow[i - 1]:
                run_lwbmatch_flow[i] = True

        [self.run_dia_umpire, self.run_ms_convert, self.run_db_generate, self.run_comet, self.run_xinteract,
         self.run_extract_ground_truth] = run_dia_umpire_flow

        [self.run_feature_finder, self.run_lwbmatch] = run_lwbmatch_flow

        if self.run_lwbmatch:
            self.run_plot_figure = True

        if run_dia_umpire_flow[-1] or run_lwbmatch_flow[-1]:
            self.run_recall_precision = True


if __name__ == '__main__':
    wfc = WorkflowController()
    wfc.run_dia_umpire = True;
    wfc.reset_workflow()
    print wfc.__dict__

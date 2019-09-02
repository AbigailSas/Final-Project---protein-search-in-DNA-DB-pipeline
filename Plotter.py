import sys
import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from collections import defaultdict

IRRELEVANT_LINES = 15
INDEX_OF_SEQ_DIRECTION = 43
MONTH_ABBREVIATION = 3
EVAL_LOC = 4
SEQ_NAME_LOC = 1
SCORE_LOC = 2
MONTH_TO_NUMBER_DICT = {"jan":1,"feb":2,"mar":3,"apr":4,"may":5,"jun":6,"jul":7,"aug":8,
                                       "sep":9,"oct":10,"nov":11,"dec":12}

class Plotter:

    def amount_chl_by_month(self, normalized_f,normalized_a,year ):
        """
        Plot amount of normalized chlorophyll along 1 year, with y axis as month and x axis as amount. 
        :param normalized_f:
        :param normalized_a:
        :param year: data from year 1 or 2
        :return: show figuers
        """
        plt.plot(normalized_f, label="chlF synthase", linewidth=2)
        plt.plot(normalized_a, label="A synthase", linewidth=2)
        plt.title(" year " + str(year), fontsize=20)
        plt.xlabel("months", fontsize=16)
        plt.ylabel("normalized chlF synthase amount", fontsize=16)
        plt.legend()
        plt.xticks(np.arange(12) + 0.25 / 2,
                   ("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"))
        plt.show()



    def relation_between_chlf_and_chla(self, f, a, year):
        """
        plot relation between amount of chlF and chlA
        :param f:
        :param a:
        :param year: data from year 1 or year 2
        :return:
        """
        months = MONTH_TO_NUMBER_DICT.keys()
        y_pos = np.arange(len(months))
        f_vals = f
        a_vals = a

        if year == 2:
            f_vals[8] = 1
            a_vals[8] = 1
        mpl_fig = plt.figure()
        percent_plot = mpl_fig.add_subplot(111)

        totals = [i + j for i, j in zip(f_vals, a_vals)]
        f_percent = [i / j * 100 for i, j in zip(f_vals, totals)]
        a_percent = [i / j * 100 for i, j in zip(a_vals, totals)]

        percent_plot.bar(range(len(f_percent)), f_percent, color='forestgreen', label="ChlF")
        percent_plot.bar(range(len(a_percent)), a_percent, bottom=f_percent, color='darkgreen', label="ChlA")

        plt.sca(percent_plot)
        plt.xticks(y_pos, months)
        percent_plot.set_ylabel("Percentage")
        percent_plot.set_xlabel("month")
        plt.title("relation between amount of chlF and chlA")
        plt.legend()
        plt.show()


    def correlation_f_a_both_years(self, normalized_f_year_1, normalized_a_year_1, normalized_f_year_2,
                                   normalized_a_year_2):
        """

        :param normalized_f_year_1:
        :param normalized_a_year_1:
        :param normalized_f_year_2:
        :param normalized_a_year_2:
        :return:
        """
        plt.scatter(normalized_f_year_1, normalized_a_year_1, label="year 1")
        plt.scatter(normalized_f_year_2, normalized_a_year_2, label="year 2")
        plt.xlabel("chl F sequences")
        plt.ylabel("chl A sequences")
        plt.legend()
        plt.show()

import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

class Visualization(object):

    def __init__(self, command, data, outfile):
        self.blue_palette = "PRGn"
        self.green_palette = "BuGn_d"
        self.data = data
        self.command = command
        self.outfile = outfile

        self._start()

    def _start(self):
        if command == "barplot":
            self.freq_barplot()
        elif command == "lmplot":
            self.lmplot()
        elif command == "disease":
            self.lmplot_by_disease()
        elif command == "viral":
            self.viral_load()

    def viral_load(self):
        sns.set(style="darkgrid")
        x = [[1, 2, 3], [4, 5], [6, 7, 8], [9, 10, 11],
             [12, 13, 14], [15, 16, 17, 18], [19, 20, 21], [22, 23, 24]]
        y = [[380000, 123000, 5110], [83000, 3360], [322000, 11200, 91.1], [98200, 3520, 47.2],
             [807000, 17600, 1520], [232000, 19300, 2230, 94], [40500, 26100, 111],
             [257000, 39200, 5850]]
        for i in range(0, len(x)):
            data = pd.DataFrame(np.array([[x1, y1] for x1, y1 in zip(x[i], y[i])]),
                                columns=["sample", "viral"])
            #plt.plot("sample", "viral", data,
            #         palette=self.green_palette,
            #         marker="-")
            plt.plot(x[i], y[i], "o-", linewidth=3, ms=12)
            plt.xlabel("Sample")
            plt.ylabel("Viral Load")
            plt.title("Viral Load of Patient " + str(i+1))

            plt.savefig("/data/chenga/regplot" + str(i+1) + ".png", dpi=200)
            plt.clf()

    def freq_barplot(self):
        sns.set(style="darkgrid")

        #x = np.array([1, 2, 3])
        #y = []
        #y.append(np.array([5.51, 5.91, 5.25]))
        #y.append(np.array([5.66, 5.03, 5.57]))
        #y.append(np.array([4.08, 4.39, 3.98]))
        #y.append(np.array([7.22, 6.15, 4.25]))

        for i in range(1, 9):
            pos = self.data[self.data.patient == i].groupby("position").groups
            sorted_pos = sorted(pos.items(), key=lambda q: q[0])

            x = np.array([int(self.data.ix[v].sample) for v in sorted_pos[0][1]])
            y = []
            z = []
            for k1, v1 in sorted_pos:
                z.append(k1)
                y.append(np.array([self.data.ix[v2].hapfreq_p for v2 in v1]))

            f, ax = plt.subplots(len(z), 1, figsize=(4,6), sharex=True)
            for n, (axis, pos) in enumerate(zip(ax, z)):
                sns.barplot(x, y[n], palette=self.green_palette, ax=axis)
                axis.plot(range(len(x)), y[n], "ko-", ms=5, mew=3)
                for label in axis.get_xticklabels() + axis.get_yticklabels():
                    label.set_fontsize(8)

                axis.set_ylim(0, max(y[n])+5.0)

                ax[n].set_ylabel("Pos " + str(pos), fontsize=8)

            f.suptitle("Patient " + str(i) + " - " \
                       + "Timeseries barplots of alt. haplotype" \
                       + "\n" + "with maf greater than 5%",
                       fontsize=10)

            ax[len(z)-1].set_xlabel("Samples", fontsize=8)
            plt.savefig("/data/chenga/ts" + str(i) + ".png", dpi=200)

    def violinplot(self):
        sns.set(style="ticks")
        e = np.array([[6.1506, 5.758, 4.8478, 4.2768, 1.70639],
                      [4.2588, 3.7663, 4.1674, 3.5528, 3.4721, 6.9316, 5.2319, 0.2309],
                      [2.0883, 4.3466, 2.5628, 6.0265, 5.6181, 0.00]])
        e.offset_spines()
        sns.violinplot(e)
        sns.despine(trim=True)
        plt.savefig(out, dpi=200)

    def boxplot(self):
        sns.set(style="darkgrid")
        data = pd.read_csv(filename)
        g = sns.factorplot("day", "freq_percent", data=data,
                           kind="box", palette=self.blue_palette, aspect=1.25)
        g.despine(left=True)
        g.set_axis_labels("Day", "Haplotype Frequency")
        plt.savefig(out, dpi=200)

    def lmplot(self):
        sns.set(style="darkgrid")
        g = sns.lmplot("day", "hapfreq_p", data,
                       hue="pos",
                       col="pos",
                       aspect=1.25,
                       hue_order=[2810, 2719],
                       legend=False,
                       legend_out=False,
                       palette="deep")
        g.set_xlabels("Days to Effervescence")
        g.set_ylabels("Alt. Haplotype Frequency (%)")
        #g.set_titles("Linear regression plot of alternate haplotypes")

        plt.savefig(self.outfile, dpi=200)

    def lmplot_by_disease(self):
        sns.set(style="darkgrid")
        data = self.data[self.data["pos"] == 2810]
        g = sns.lmplot("day", "hapfreq_p", data,
                       hue="disease",
                       aspect=1.25,
                       palette="deep")
        g.set_xlabels("Days to Effervescence")
        g.set_ylabels("Alt. Haplotype Frequency (%)")

        plt.savefig(self.outfile, dpi=200)

if __name__ == "__main__":

    command, csv_file, outfile = sys.argv[1:]
    
    data = pd.read_csv(csv_file)
    v = Visualization(command, data, outfile)

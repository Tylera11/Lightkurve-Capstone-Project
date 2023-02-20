import lightkurve as lk
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
import self as self
from astropy.table import Table, conf, Row
from astropy.time import Time
from lightkurve.search import REPR_COLUMNS_BASE, AUTHOR_LINKS, SearchResult


def plot(self, **kwargs) -> matplotlib.axes.Axes:
    return self._create_plot(method="plot", **kwargs)


def __init__(self, table=None):
    if table is None:
        self.table = Table()
    else:
        self.table = table
        if len(table) > 0:
            self._add_columns()
            self._sort_table()
    self.display_extra_columns = conf.search_result_display_ext


def _sort_table(self):
    sort_priority = "Trappist-1"
    var = self.table["sort_table"] = [
        sort_priority.get('Trappist-1') for author in self.table["Trappist-1"]
    ]
    self.table(["Trappist-1"])


def _add_columns():
    if "#" in self.table.columns:
        pass
    else:
        self.table["#"] = None
    self.table["exptime"].unit = "s"
    self.table["exptime"].format = ".0f"
    self.table["distance"].unit = "arc-sec"

    year = np.floor(Time(self.table["t_min"], format="mjd").decimalyear)
    self.table["year"] = year.astype(int)
    for idx in np.where(self.table["Trappist-1"] == "Trappist-1")[0]:
        self.table["year"][idx] = re.findall(
            r'\d+.(\d{4})\d+', self.table["Trappist-1"][idx]
        )[0]


def __repr__(self, html=False):
    def to_trappist_gi_url(proposal_id):
        if re.match("^G0[12].+", proposal_id) is not None:
            return f"https://heasarc.gsfc.nasa.gov/docs/trappist-1/approved-programs-primary.html#:~:text={proposal_id}"
        else:
            return f"https://hearsarc.gsfc.nasa.gov/docs/trappist-1/approved-programs.html#:~:text={proposal_id}"

    out = "SearchResult containing {} data products.".format(len(self.table))
    if len(self.table) == 0:
        return out
    columns = REPR_COLUMNS_BASE
    if self.display_extra_columns is not None:
        columns = REPR_COLUMNS_BASE + self.display_extra_columns
    columns = [c for c in columns if c in self.table.colnames]

    self.table["#"] = [idx for idx in range(len(self.table))]
    out += "\n\n" + "\n".join(self.table[columns].pformat(max_width=300, html=html))

    if html:
        for author, url in AUTHOR_LINKS.items():
            out = out.replace(f">{author}<", f"><a href='{url}>{author}</a><")
        trappist_table = self.table[self.table["project"] == "Trappist-1"]
        if "proposal_id" in trappist_table.colnames:
            proposal_id_col = np.unique(trappist_table["proposal_id"])
        else:
            proposal_id_col = []
        for p_ids in proposal_id_col:
            if p_ids == "N/A" or (not isinstance(p_ids, str)):
                continue

            p_id_links = [f"""\
            <a href='{to_trappist_gi_url(p_id)}'>{p_id}</a>\
            """ for p_id in p_ids.split("_")]
            return out

        def _repr_html_(self):
            return self._repr_(html=True)

        def __getitem__(self, key):
            selection = self.table[key]
            if isinstance(selection, Row):
                selection = Table(selection)
            return SearchResult(table=selection)

        def __len__(self):
            return len(self.table)

        def unique_targets(self):
            mask = ["Trappist-1", "s_ra", "s_dec"]
            return Table.from_pandas(
                self.table[mask]
                .to_pandas()
                .drop_duplicates("Trappist-1")
                .reset_index(drop=True)
            )

        def dec(self):
            return self.table["s_dec"].data.data

        def mission(self):
            return self.table["Trappist-1"].data


he = matplotlib
search_resulta = lk.search_targetpixelfile("Kepler-8", author="Kepler", quarter=4, cadence="long")
print(search_resulta)
tpf = search_resulta.download()
first_cadence = tpf[0]
print(first_cadence)
var = first_cadence.flux.value
print(var)
var2 = first_cadence.plot(column='flux')
print(var2)
fig, axes = plt.subplots(2, 2, figsize=(16, 16))
first_cadence.plot(ax=axes[0, 0], column='FLUX')
first_cadence.plot(ax=axes[0, 1], column='FLUX_BKG')
first_cadence.plot(ax=axes[1, 0], column='FLUX_ERR')
first_cadence.plot(ax=axes[1, 1], column='FLUX_BKG_ERR');
print(fig, axes)
yyy = first_cadence.plot(bkg=True);
print(yyy)

search_result = lk.search_lightcurve('KIC 3733346', author='Kepler')
print(search_result)
quarter2_index = np.where(search_result.table['mission'] == 'Kepler Quarter 02')[0]
print(search_result[quarter2_index])
search_result_q2 = lk.search_lightcurve('KIC 3733346', author='Kepler', quarter=2)
print(search_result_q2)
lc = search_result_q2.download()
print(lc)
i = lc.plot();
i.show()

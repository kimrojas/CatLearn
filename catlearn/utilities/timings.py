from time import time
from datetime import datetime

import sys
import numpy as np
import tabulate

tabulate.PRESERVE_WHITESPACE = True
from tabulate import tabulate


class Timer:
    def __init__(
        self,
        debug=False,
        unit="s",
        floatfmt=".6f",
        logfile=sys.stdout,
    ):
        self.debug = debug
        self.unit = unit.lower()
        self.floatfmt = floatfmt
        self.logfile = logfile

        self.timings = {}
        self.unit_mult = {"s": 1, "m": 1 / 60, "h": 1 / 3600}[self.unit]
        self.tab = "    "

    def timelog(self, msg, mode="w"):
        if self.logfile == sys.stdout:
            print(msg)
        if self.logfile is None:
            return

        with open(self.logfile, mode) as f:
            f.write(msg + "\n")

    @staticmethod
    def get_time():
        _datetime = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        _time = time()
        return _datetime, _time

    def start(self, name, level=0):
        assert name not in self.timings, "Timer already started"

        _datetime, _time = self.get_time()
        self.timings[name] = {
            "name": name,
            "start_time": _time,
            "start_datetime": _datetime,
            "end_time": None,
            "end_datetime": None,
            "duration": None,
            "level": level,
        }

        if self.debug:
            self.timelog(f"[ {_datetime} ] START  | {name} ...")

    def stop(self, name):
        assert name in self.timings, "Timer not started"

        _datetime, _time = self.get_time()
        t_obj = self.timings[name]
        t_obj["end_time"] = _time
        t_obj["end_datetime"] = _datetime
        t_obj["duration"] = (t_obj["end_time"] - t_obj["start_time"]) * self.unit_mult

        if self.debug:
            self.timelog(
                f"[ {_datetime} ] STOP   | {name} ... {t_obj['duration']:.6f} {self.unit}"
            )

        self.report()

    def report(self, ml_iter=False):
        headers = ["Task Name", "Start Time", "End Time", f"Duration ({self.unit})"]
        table = []
        for k, v in self.timings.items():
            if not ml_iter and "ML NEB - cycle" in k:
                continue
            taskname = self.tab * v["level"] + v["name"]
            row = [
                taskname,
                v["start_datetime"],
                v["end_datetime"],
                v["duration"],
            ]
            table.append(row)
        kwargs = dict(
            floatfmt=".6f",
            tablefmt="grid",
            numalign="decimal",
        )

        _table = tabulate(table, headers, **kwargs)

        self.timelog(_table)

        return _table

    def summarize(self, patterns):
        total_time = self.timings["MLNEB Run"]["duration"]

        # ----  Summary logging ----
        def summary_log():
            data = {p: {"name": p, "data_list": []} for p in patterns}

            for k, v in self.timings.items():
                for p in patterns:
                    if p in k:
                        data[p]["data_list"].append(v["duration"])

            for k, v in data.items():
                v["mean"] = np.mean(v["data_list"])
                v["std"] = np.std(v["data_list"])
                v["min"] = np.min(v["data_list"])
                v["max"] = np.max(v["data_list"])
                v["count"] = len(v["data_list"])
                v["sum"] = np.sum(v["data_list"])
                v["percent"] = v["sum"] / total_time * 100

            headers = [
                "Task Name",
                "Mean",
                "Std",
                "Min",
                "Max",
                "Count",
                "Runtime",
                "% Runtime Share",
            ]
            table = []
            for k, v in data.items():
                row = [
                    v["name"],
                    v["mean"],
                    v["std"],
                    v["min"],
                    v["max"],
                    v["count"],
                    v["sum"],
                    f"{v['percent']:.2f}",
                ]
                table.append(row)

            kwargs = dict(
                floatfmt=(None, ".6f", ".6f", ".6f", ".6f", "d", ".6f", ".2f"),
                tablefmt="grid",
                numalign="decimal",
            )

            _table = tabulate(table, headers, **kwargs)

            self.timelog(_table, mode="a")

            return data

        # ----  MDMin iter logging ----

        def MDMin_iter_log():
            pattern = "ML NEB - cycle"
            patterndict = {k: v for k, v in self.timings.items() if pattern in k}
            datadict = {}
            for k, v in patterndict.items():
                iter = int(k.split()[0].strip("()"))
                cycle = int(k.split()[-1].strip("[]"))
                if iter not in datadict:
                    datadict[iter] = {}
                    datadict[iter]["data_list"] = []
                datadict[iter]["data_list"].append(v["duration"])

            for k, v in datadict.items():
                v["mean"] = np.mean(v["data_list"])
                v["std"] = np.std(v["data_list"])
                v["min"] = np.min(v["data_list"])
                v["max"] = np.max(v["data_list"])
                v["count"] = len(v["data_list"])
                v["sum"] = np.sum(v["data_list"])
                v["percent"] = v["sum"] / total_time * 100

            headers = [
                "ML Iter",
                "Mean",
                "Std",
                "Min",
                "Max",
                "Count",
                "Runtime",
                "% Runtime Share",
            ]
            table = []
            for k, v in datadict.items():
                row = [
                    k,
                    v["mean"],
                    v["std"],
                    v["min"],
                    v["max"],
                    v["count"],
                    v["sum"],
                    v["percent"],
                ]
                table.append(row)
            kwargs = dict(
                floatfmt=(None, ".6f", ".6f", ".6f", ".6f", "d", ".6f", ".2f"),
                tablefmt="grid",
                numalign="decimal",
            )
            _table = tabulate(table, headers, **kwargs)
            self.timelog(_table, mode="a")

            return datadict

        # -- Call functions --
        self.data_summary = summary_log()
        self.data_mdmin = MDMin_iter_log()

    def timer_decorator(self, name):
        def decorator_func(func):
            def wrapper(*args, **kwargs):
                self.start(name)
                result = func(*args, **kwargs)
                self.stop(name)
                return result

            return wrapper

        return decorator_func

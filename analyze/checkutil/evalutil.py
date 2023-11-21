import pathlib
import re

# Regexes
patt1 = re.compile(r"bin/sbcs\.linux\.x86_64\.gnu\.opt\.spx2 (\d+)")
patt2 = re.compile(r"Constraint Timings :")
patt3 = re.compile(r"  (symresack|symretope|orbisack|orbitope)\s+:\s+(\d+(?:\.\d+)?)\+?\s+(\d+(?:\.\d+)?)\+?\s+(\d+(?:\.\d+)?)\+?\s+(\d+(?:\.\d+)?)\+?\s+(\d+(?:\.\d+)?)\+?\s+(\d+(?:\.\d+)?)\+?\s+(\d+(?:\.\d+)?)\+?\s+(\d+(?:\.\d+)?)\+?\s+(\d+(?:\.\d+)?)\+?\s+(\d+(?:\.\d+)?)\+?")
patt4 = re.compile(r"SCIP Status\s+: ([^\[]+)\[([^\[]+)\]")
patt5 = re.compile(r"Primal Bound\s+: ([+-]?(\d+([.]\d*)?(e[+-]?\d+)?|[.]\d+(e[+-]?\d+)?))")
patt6 = re.compile(r"Dual Bound\s+: ([+-]?(\d+([.]\d*)?(e[+-]?\d+)?|[.]\d+(e[+-]?\d+)?))")
patt7 = re.compile(r"@01\s*([^=\s]*)\s*=+$")
patt8 = re.compile(r"@03 (\d+)")
patt9 = re.compile(r"@04 (\d+)")
patt10 = re.compile(r"Solving Time \(sec\)\s+:\s+(\d+\.?\d+)")
patt11 = re.compile(r"(?:reading parameter file|reading user parameter file|Reading parameters from) <([^>]*)>(?: ...)?")


class Instance:
    def __init__(self, *args) -> None:
        (self.outpath, self.instancepath, self.settingspath, self.settings, self.instancename, self.totaltime,
            self.constraint_times, self.solved, self.status, self.primal, self.dual, self.objdir, self.crashed) = args


suffixes = [".gz", ".mps", ".cip"]
def remove_suffices(s):
    has_suff = True
    while has_suff:
        has_suff = False
        for suff in suffixes:
            if s[-len(suff):] == suff:
                s = s[:-len(suff)]
                has_suff = True
                continue
    return s


def read_outfile(outpath : pathlib.Path):
    instance = None
    totaltime = None
    timings_section = False
    timings_times = {}
    solved = None
    status = None
    primal = None
    dual = None
    objdir = None
    settingsfile = None
    instancename = None
    starttime = None
    endtime = None

    with open(outpath, "r") as f:
        for line in f:
            m = patt7.match(line)
            if m is not None:
                # Save previous instance
                if instance is not None:
                    yield Instance(outpath, instance, settingsfile, settingsfile.stem, instancename, totaltime,
                        timings_times.copy(), solved, status, primal, dual, objdir)
                    instance = None
                    totaltime = None
                    timings_section = False
                    timings_times.clear()
                    solved = None
                    status = None
                    primal = None
                    dual = None
                    objdir = None
                    starttime = None
                    endtime = None

                instance = pathlib.Path(m.group(1))
                instancename = remove_suffices(instance.name)

                totaltime = None
            m = patt8.match(line)
            if m is not None:
                starttime = int(m.group(1))
            m = patt9.match(line)
            if m is not None:
                endtime = int(m.group(1))

            m = patt10.match(line)
            if m is not None:
                totaltime = float(m.group(1))
            m = patt11.match(line)
            if m is not None:
                settingsfile = pathlib.Path(m.group(1))

            # Hamilton situation!
            if instance == pathlib.Path(""):
                m = patt1.match(line)
                if m is not None:
                    instance = pathlib.Path(f"hamilton{m.group(1)}")
                    instancename = m.group(1)

            m = patt2.match(line)
            if m is not None:
                timings_section = True
            if timings_section:
                m = patt3.match(line)
                if m is not None:
                    constraint = m.group(1)
                    if constraint not in timings_times:
                        # print(line)
                        ctrtotaltime = float(m.group(2))
                        proptime = float(m.group(5))
                        resproptime = float(m.group(10))
                        timings_times[constraint] = (ctrtotaltime, proptime, resproptime)
                        # print(constraint, ctrtotaltime, proptime, resproptime)

            m = patt4.match(line)
            if m is not None:
                solved = m.group(1).strip()
                status = m.group(2)
            m = patt5.match(line)
            if m is not None:
                primal = float(m.group(1))
            m = patt6.match(line)
            if m is not None:
                dual = float(m.group(1))
            # if objdir is None:
            #     # Take the first one, which is the objective direction of the actual problem (not the presolved one)
            #     m = re.match(r"  Objective\s+: (\w+),", line)
            #     if m is not None:
            #         objdir = m.group(1)

        # Deal with instances exceeding time limit.
        if totaltime is None:
            if totaltime is None and starttime is not None and endtime is not None and endtime - starttime > 7200:
                totaltime = 7200.0
                solved = "stuck"
                status = "time limit reached"
                primal = -1.0
                dual = -1.0
                objdir = None

        if instance is not None:
            crashed = any(i is None for i in [outpath, instance, settingsfile, instancename, totaltime, solved, status, primal, dual])
            if totaltime is None:
                totaltime = 7200 # Time limit
            if solved is None:
                assert status is None
                solved = "stuck"
                status = "time limit reached"

            yield Instance(outpath, instance, settingsfile, settingsfile.stem, instancename, totaltime, 
                timings_times.copy(), solved, status, primal, dual, objdir, crashed)
            instance = None
            totaltime = None
            timings_section = False
            timings_times.clear()
            solved = None
            status = None
            primal = None
            dual = None
            objdir = None
            starttime = None
            endtime = None


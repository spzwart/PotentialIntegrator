from amuse.units import units

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", 
                      dest="filename", 
                      default = "binary.amuse",
                      help="initial conditions datafile [%default]")
    result.add_option("--fig", 
                      dest="figname", 
                      default = "orbit.pdf",
                      help="figure filename [%default]")
    result.add_option("-t", type="float", unit=units.yr,
                      dest="t_end", 
                      default = 10|units.yr,
                      help="end time of the simulation [%default]")
    result.add_option("--dt", type="float", unit=units.yr,
                      dest="dt", 
                      default = 0.25|units.yr,
                      help="time stel [%default]")
    return result

if __name__ == "__main__":
    o, arguments = new_option_parser().parse_args()
    print(o)
    print(arguments)

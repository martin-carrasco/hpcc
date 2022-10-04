import pandas as pd


HEADERS = ["Type", "Rank", "Container", "Start", "End", "Duration", "Level", "State"]
def parse():
    df = pd.read_csv('smpi_simgrid.trace', names=HEADERS)
    print(df)

if __name__ == '__main__':
    parse()
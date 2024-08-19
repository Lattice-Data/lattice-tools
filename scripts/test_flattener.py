import lattice
import multiprocessing
import traceback
from flattener import main


mode = "demo"
connection = None

def set_global_connection():
    global connection
    if not connection:
        connection = lattice.Connection(mode)


def flatten(file):
    try:
        main(file, connection)
        return (file, "SUCCESS")
    except Exception as e:
        return (file, e)


files = [
    "LATDF190KNY",
    "LATDF584NGT",
    "LATDF742BQI",  # np.bool error fixed on main, not here yet
    "LATDF329OGL",
    "LATDF477OUM",
    "LATDF483VSD",
    "LATDF994BQY",
    "LATDF448PLU",
    "LATDF366MLP",
    "LATDF216UIK",  # fix merged on main, not yet on this branch
]

if __name__ == "__main__":
    results = list()
    with multiprocessing.Pool(initializer=set_global_connection) as pool:
        iterator = pool.imap(flatten, files)
        while True:
            try:
                results.append(next(iterator))
            except StopIteration:
                break
            except Exception as e:
                print(e)
                break

    print("FINAL RESULTS:")
    print("=" * 80)
    for file, returned_obj in results:
        if isinstance(returned_obj, str):
            print(f"{file}: SUCCESS")
        else:
            print(f"{file}: FAILURE")
            traceback.print_exception(None, returned_obj, returned_obj.__traceback__)
            print("=" * 80)

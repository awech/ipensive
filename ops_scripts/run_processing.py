import time
import logging
from obspy import UTCDateTime as utc
from ipensive import ipensive_utils as utils
from ipensive import array_processing


def run_array_processing():
    """
    Main entry point for the script. Handles argument parsing, configuration loading,
    and processing of arrays.
    """
    timer_0 = time.time()

    args = array_processing.parse_args()  # Parse command-line arguments
    config_file = args.config
    config = utils.load_config(config_file)  # Load configuration

    config["plot"] = False if args.no_plot else True

    T0, delay = array_processing.get_starttime(config, args)  # Determine start time and delay

    utils.setup_logging(T0, config, arg_opt=args.log)
    my_log = logging.getLogger(__name__)
    my_log.info("Process Started")
    
    time.sleep(delay)  # Pause to allow for data latency to catch up
    timer_0 += delay

    if args.array:
        # Process a single array if specified
        timer_tmp = time.time()
        array_processing.process_array(config, args.array.replace("_", " "), T0)
        dt = time.time() - timer_tmp
        my_log.info(f"{dt:.1f} seconds to process {args.array.replace('_', ' ')}")
    else:
        # Process all arrays in the configuration
        for array_name in config["array_list"]:
            timer_tmp = time.time()
            array_processing.process_array(config, array_name, T0)
            dt = time.time() - timer_tmp
            my_log.info(f"{dt:.1f} seconds to process {array_name}\n\n")

    # Write out the new HTML file
    array_processing.write_html(config)

    my_log.info(f"{time.time()-timer_0:.1f} seconds to process all")

if __name__ == "__main__":
    run_array_processing()

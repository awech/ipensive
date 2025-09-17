import time
import logging
from ipensive import ipensive_utils as utils
from ipensive import array_processing


if __name__ == "__main__":
    """
    Main entry point for the script.
    """
    
    timer_0 = time.time()

    args = array_processing.parse_args()  # Parse command-line arguments
    ipensive_config_file = args.config
    config = utils.load_ipensive_config(ipensive_config_file)  # Load configuration

    config["plot"] = False if args.no_plot else True

    T0, delay = array_processing.get_starttime(config, args)  # Determine start time and delay

    utils.setup_logging(T0, config, arg_opt=args.log)
    my_log = logging.getLogger(__name__)
    my_log.info(f"Using ipensive config file: {ipensive_config_file}")
    for f in config["array_config_files"]:
        my_log.info(f"Using array config file: {f}")
    my_log.info("Process Started")

    if delay > 0:
        my_log.info(f"Waiting {delay:.1f} seconds for data latency to catch up")
        time.sleep(delay)
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
            try:
                array_processing.process_array(config, array_name, T0)
            except Exception as ex:
                my_log.error(f"Error processing {array_name}: {ex}")
            dt = time.time() - timer_tmp
            my_log.info(f"{dt:.1f} seconds to process {array_name}\n")

    # Write out the new HTML file
    utils.write_html(config)

    my_log.info(f"{time.time() - timer_0:.1f} seconds to process all")
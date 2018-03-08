""" license 
"""

from __future__ import print_function # use print as a function

import sys

from utils import logger # logger facilitates acces to debug, warning, error...

try:
    if hasattr(sys, "run_from_cmdl") is True: #hasattr function return True if string is an atribute of the first object
        raise ImportError
    from pycompss.api.parameter import FILE_IN, FILE_OUT
    from pycompss.api.task import task
    from pycompss.api.api import compss_wait_on
except ImportError:
    logger.warn("[Warning] Cannot import \"pycompss\" API packages.")
    logger.warn("Using mock decorators")

    from utils.dummy_pycompss import FILE_IN, FILE_OUT # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import task # pylint: disable=ungrouped-imports
    from utils.dummy_pycompss import compss_wait_on # pylint: disable=ungrouped-imports
    

from basic_modules.tool import Tool 
"""
class basic_modules.tool.Tool(configuration = {})
    Tool is a specific operation done in a input file generating an output file
    A Tool is executed by runnign the run() method which support multiples inputs an outputs
    run() method recieved metadata for each of the input data elements. Tools responsability
    to generate the metadata for each of the output data elements, returned in a tuple.

    run method calls the methods that perform the operations to implements the Tool's functionality
    each method should be decorated using the @task decorator. Task constraints can be configured 
    using the @constraint decorator.

    run(input_files, input_metadata, output_files) 
        perform the functionality of the Tool. Usually > 
        1)Import specific tool libraries,
        2) check the input data, 
        3) convert input data into tool-specific format, 
        4) Perform the tool specific operations, 
        5) convert output into output format, 
        6) write metadata for output,
        7) handle error for any of the above

        If failure return None with the exception instance to the 
        output metadata. 

        The method will run the task, each task shoudl have unique name 
        that identify the operation. used then by COMPSs runtime to build
        a graph and trace.

        parameters
            input_file (dict), a dict of absolute path names of the input data elements, associ- ated with their role;
            input_metadata (dict), a dict of metadatas for each of the input data elements, associated with their role;
            output_files (dict), a dict of absolute path names of the output data elements, associated with their role.

        Returns
            output_files [dict] a dict of absolute path names of the output data elements created by the
            Tool, associated with their role;
            output_metadata [dict] a dict of metadatas for each of the output data elements created by the Tool,
            associated with their role;

        Return type (output_files, output_metadata)



"""

from basic_modules.metadata import Metadata # contain the metadata of the files 


class testToolPablo(Tool):
    """ tool for printing in a file if the number of 
        characters is even or odd
    """
    def __init__(self, configuration = None):
        """ init function"""
        logger.info("Test writer") # looger.info(" ") logs a message with a level INFO. args are interpreted as for debug()
        """ utils.logger.debug(message, *args, **kargs) 
            message is the message format string, and args are the arguments which are merged into 
            the msg using the string formatting operator.
        """
        Tool.__init__(self)

        if configuration is None:
            configuration = {}

        self.configuration.update(configuration)


    @task(returns=bool, file_in_loc=FILE_IN, file_out_loc=FILE_OUT, isModifier=False)
    def writterEvenOrOdd(self, file_in_loc, file_out_loc):
        """
        Count character on a file and return even or odd

        ----------------------------------------------------

        Parameters
        ----------

            file_in_loc: str
                location input file
            file_out_loc:str
                location output file

        Returns
        -------

        bool


        Writes to the file defined by COMPSs on the location
        """

        try:
            with open(file_in_loc, "r") as file:
                with open(file_out_loc, "w") as handle_out:
                    total_char = 0
                    for line in file:
                        total_char += len(line)

                    if total_char%2 == 0:
                        handle_out.write("The number of characters is even")
                    else:
                        handle_out.write("The number of characters is odd")

        except IOError as error:
            logger.fatal("I/O error({0}): {1}".format(error.errno, error.strerror))
            return False

        return True

    def run(self, input_files, input_metadata, output_files):
        """
        Run function that run the function writtenEvenOrOdd

        Parameters:
        -----------

        input_file: dict
            list of input files - in case there is non

        input_metadata: dict 
            metadata matching the input files

        output_files: dict
            list output files

        Returns
        -------

        output_files: dict

        output_metadata: dict
            metadata that match with the input files
        """

        results = self.writtenEvenOrOdd(input_files["inpput"],
                                        output_files["output"])

        results = compss_wait_on(results)

        if results == False:
            logger.fatal("Test writterEvenOrOdd: caca de la vaca paca")
            return {}, {}

        output_metadata = {
            "output": Metadata(
                data_type="<data_type>",
                file_type="txt",
                file_path=output_files["output"],
                sources=[input_metadata["input"].file_path],
                taxon_id=input_metadata["input"].taxon_id,
                meta_data={
                    "tool": "testTool"
                }
            )
        }
        
        return (results, output_metadata)

test1 = testToolPablo()

characters_result = test1.writterEvenOrOdd("/Users/pacera/test_pipeline/mg-process-test1/tool/data_test.txt","/Users/pacera/test_pipeline/mg-process-test1/tool/output_file.txt")

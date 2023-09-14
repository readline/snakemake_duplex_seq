

def ln(files, outdir):
    """Creates symlinks for files to an output directory.
    @param files list[<str>]:
        List of filenames
    @param outdir <str>:
        Destination or output directory to create symlinks
    """
    # Create symlinks for each file in the output directory
    for file in files:
        ln = os.path.join(outdir, os.path.basename(file))
        if not exists(ln):
                os.symlink(os.path.abspath(os.path.realpath(file)), ln)
                
                
def safe_copy(source, target, resources = []):
    """Private function: Given a list paths it will recursively copy each to the
    target location. If a target path already exists, it will NOT over-write the
    existing paths data.
    @param resources <list[str]>:
        List of paths to copy over to target location
    @params source <str>:
        Add a prefix PATH to each resource
    @param target <str>:
        Target path to copy templates and required resources
    """

    for resource in resources:
        destination = os.path.join(target, resource)
        if not exists(destination):
            # Required resources do not exist
            copytree(os.path.join(source, resource), destination)
            
            
def provided(samplelist, condition):
    """
    Determines if optional rules should run. If an empty list is provided to rule all,
    snakemake will not try to generate that set of target files. If a given condition
    is not met (i.e. False) then it will not try to run that rule.
    """

    if not condition:
        # If condition is False, 
        # returns an empty list 
        # to prevent rule from 
        # running
        samplelist = []

    return samplelist


def ignore(samplelist, condition):
    """
    Determines if optional rules should run. If an empty list is provided to rule all,
    snakemake will not try to generate that set of target files. If a given condition
    is met (i.e. True) then it will not try to run that rule. This function is the 
    inverse to provided(). 
    """

    if condition:
        # If condition is True, 
        # returns an empty list 
        # to prevent rule from 
        # running
        samplelist = []

    return samplelist
import ConfigParser

def Config(config_file,section):
    config = ConfigParser.ConfigParser()
    try:
        config.readfp(open(config_file))
    except:
        message = "could not find the config file " + config_file
        exit(message)

    Version = config.get(section,'Version')
    Location = config.get(section,'Location')
    return(Version,Location)

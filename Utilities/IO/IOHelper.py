import os
from logging import error


abs_file = os.path.abspath(os.path.dirname(__file__))
MAIN_DIR = os.sep.join(abs_file.split(os.sep)[:-2]) + os.sep
CONFIG_FILE = "config.txt"
CONFIGT_FILE = "configt.txt"
CAMCONFIG_FILE = "cameraconfig.txt"
CONFIG_FILE_PATH = MAIN_DIR + CONFIG_FILE
CONFIGT_FILE_PATH = MAIN_DIR + CONFIGT_FILE
CAMCONFIGT_FILE_PATH = MAIN_DIR + CAMCONFIG_FILE


DATA_PATH = "DATA_PATH"

#modi creatconfig
def create_config_file(config_path=CONFIG_FILE_PATH):
    """
    generate a generic config file from a given path
    :param config_path:
    :return:
    """
    of = open(config_path, 'w')
    # of.write("{}={}Debug\images{}\n".format(DATA_PATH, MAIN_DIR, os.sep))
    # AtomImgSCNU\Model\Instruments\Camera\TUCamdll\pic\ODimage
    of.write("{}={}Model\Instruments\Camera\TUCamdll\pic\ODimage{}\n".format(DATA_PATH, MAIN_DIR, os.sep))
    of.close()

def create_configt_file(config_path=CONFIGT_FILE_PATH):
    """
    generate a generic config file from a given path
    :param config_path:
    :return:
    """
    of = open(config_path, 'w')
    of.write("{}={}images\n".format(DATA_PATH, MAIN_DIR))
    of.close()

def create_camconfig_file(config_path=CAMCONFIGT_FILE_PATH):
    """
    generate a generic config file from a given path
    :param config_path:
    :return:
    """
    if os.path.exists(config_path):
        pass
    else:
        of = open(config_path, 'w')
        # of.write("{}={}images\n".format(DATA_PATH, MAIN_DIR))
        of.close()

def set_config_setting(setting, setting_value, config_file_path=CONFIG_FILE_PATH):
    """
        sets a setting to a given value inside the configuration file
    """
    try:
        # open the file
        config_file = open(config_file_path, 'r')
        # read the lines into a list
        lines = config_file.readlines()
        # loop through all the lines
        for i, line in enumerate(lines):
            # the setting should have a = sign in the line
            try:
                [left, right] = line.split("=")
                # separate the setting name from its value
                left = left.strip()  # name
                right = right.strip()  # value
                # if the name corresponds to the setting we want, we write the value
                if left == setting:
                    lines[i] = "{}={}\n".format(setting, setting_value)
            except ValueError as e:
                if "need more than 1 value to unpack" in e:
                    pass
                else:
                    raise e
        config_file.close()
        # reopen the file and write the modified lines
        config_file = open(config_file_path, 'w')
        config_file.writelines(lines)
        config_file.close()
        print(("The parameter {} in the config file was successfully changed to {}".format(setting, setting_value)))
    except:
        error("Could not set the parameter {} to {} in the config file located at {}\n".format(setting, setting_value, config_file_path))


def get_config_setting(setting, config_file_path=CONFIG_FILE_PATH):
    """
    get a setting from the configuration file
    :param setting:
    :param config_file_path:
    :return:
    """
    value = None
    try:
        config_file = open(config_file_path, 'r')
        for line in config_file:
            # '#' is used as a comment flag
            if line[0] != '#':
                [left, right] = line.split('=')
                # separate the setting name from its value
                left = left.strip()  # name
                right = right.strip()  # value
                # choose the name we want, then read the value
                if left == setting:
                    value = right
        if not value:
            print("Configuration file doesn't contain {}".format(setting))
        config_file.close()
    except IOError:
        print("No configuration file {} found".format(config_file_path))
        value = None
    return value

def get_configt_setting(setting, config_file_path=CONFIGT_FILE_PATH):
    """
    get a setting from the configuration file
    :param setting:
    :param config_file_path:
    :return:
    """
    value = None
    try:
        config_file = open(config_file_path, 'r')
        for line in config_file:
            # '#' is used as a comment flag
            if line[0] != '#':
                [left, right] = line.split('=')
                # separate the setting name from its value
                left = left.strip()  # name #Remove before and after Spaces
                right = right.strip()  # value #Remove before and after Spaces
                # choose the name we want, then read the value
                if left == setting:
                    value = right
        if not value:
            print("Configuration file doesn't contain {}".format(setting))
        config_file.close()
    except IOError:
        print("No configuration file {} found".format(config_file_path))
        value = None
    return value

def get_camconfig_setting(config_file_path=CAMCONFIGT_FILE_PATH):
    """
    get a setting from the configuration file
    :param setting:
    :param config_file_path:
    :return:
    """
    value = []
    try:
        config_file = open(config_file_path, 'r')
        # file.readlines()
        for line in config_file.readlines():
            # '#' is used as a comment flag
            if line[0] != '#':
                [name, nickname, con] = line.split(',')
                name = name.strip()  # name #Remove before and after Spaces
                nickname = nickname.strip()  # value #Remove before and after Spaces
                con = con.strip()  # value #Remove before and after Spaces
                # choose the name we want, then read the value
                value.append(name)
                value.append(nickname)
                value.append(con)
        config_file.close()
        # print(value)
    except IOError:
        print("No configuration file {} found".format(config_file_path))
        value = None
    return value

if __name__ == '__main__':
    create_config_file(CONFIG_FILE_PATH)
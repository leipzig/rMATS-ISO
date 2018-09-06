import sys, os
import time


def format_time(fp, header, str):
    fp.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str)

def out_format_time(header, str):
    sys.stdout.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str)


def err_format_time(header, str):
    sys.stderr.write('==' + time.strftime(" %H:%M:%S-%b-%d-%Y ", time.localtime()) + '== [' + header + '] ' + str)

def fatal_format_time(header, str):
    err_format_time(header, str)
    sys.exit(1)

def exec_cmd(fp, header, cmd):
    format_time(fp, header, cmd + '\n')
    return os.system(cmd)


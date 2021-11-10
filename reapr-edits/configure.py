#!/usr/bin/python

import argparse, subprocess, os

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Configure location of required software')
    parser.add_argument('--locarna-prefix', default='', help='Installation directory of LocARNA package.')
    parser.add_argument('--rnaz-prefix', default='', help='Installation directory of RNAz package.')
    parser.add_argument('--alistat', default='', help='Path to alistat command from squid package.')
    parser.add_argument('--compalignp', default='', help='Path to compalignp command.')
    args = parser.parse_args()

    this_directory = os.path.dirname(os.path.abspath(__file__))
    c = open(os.path.join(this_directory, 'commands.py'), 'w')
    c.write("mlocarna='%s'\n" % ('mlocarna' if args.locarna_prefix=='' else os.path.join(args.locarna_prefix, 'bin/mlocarna')))
    c.write("RNAz='%s'\n" % ('RNAz' if args.rnaz_prefix=='' else os.path.join(args.rnaz_prefix, 'bin/RNAz')))
    rnaz_prefix = subprocess.Popen('which RNAz', shell=True, stdout=subprocess.PIPE).communicate()[0].split('\n')[0].split('bin/RNAz')[0] if args.rnaz_prefix=='' else args.rnaz_prefix
    c.write("rnazWindow='%s'\n" % os.path.join(rnaz_prefix, 'share/RNAz/perl/rnazWindow.pl'))
    c.write("alistat='%s'\n" % ('alistat' if args.alistat=='' else args.alistat))
#    c.write("compalignp='%s'\n" % os.path.join(this_directory, 'compalignp'))
    c.write("compalignp='%s'\n" % ('compalignp' if args.compalignp=='' else args.compalignp))           
    c.close()

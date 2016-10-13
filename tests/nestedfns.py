#!/usr/bin/python2.7
# -*- coding: utf-8 -*-
"""
nestedfns.py: Test of nested functions in python

Created by Scott Feister on Thu Oct 13 15:54:55 2016

This is to be used with tstools
"""

def hello(name, lastname=None):
    """Prints 'Hello, <name>!' or 'Hello, <full name>!' """
    if not lastname:
        print("Hello, " + name + "!")
    else:
        print("Hello, " + name + " " + lastname + "!")

def bonjour(name, lastname=None):
    """Prints 'Bonjour, <name>!' or 'Bonjour, <full name>!' """
    if not lastname:
        print("Bonjour, " + name + "!")
    else:
        print("Bonjour, " + name + " " + lastname + "!")

def welcome(name, lastname=None, greeting=None):
    """ Prints a hello message (if specified) followed by a welcome message. """
    if callable(greeting): # Note: "callable()" Breaks for Python 3.0, 3.1, but works in Python 2.x, 3.2+
        greeting(name, lastname)
    print(name + ", welcome to the party!")
    
if __name__ == "__main__":
    print("First pass:")
    welcome("Bob", greeting=bonjour)
    print("Second pass:")
    welcome("Jeremy", "Fisher", greeting=hello)
    print("Third pass:")
    welcome("Stacy")

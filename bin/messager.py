"""
set of functions to print info message
"""

def showInfoMessage(message):
    print(">INFO! {}".format(message))

def showStep(message):
    print(">*** {} ***<".format(message.upper()))

def showMessageLevel(message, level=1):
    print("{0} {1}".format("*"*(4-level), message))

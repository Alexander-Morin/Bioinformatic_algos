# Chapter 1 is centered around text matching algorithms, framed by asking how to find the ori in bacteria.

def pattern_count(text, pattern):
    """Assumes text and pattern are strings. Returns the count of times a pattern was found in a string"""
    count = 0
    k = len(pattern)
    for i in range(0, len(text) - k+1):
        if text[i:i+k] == pattern:
            count =+ 1
    return count

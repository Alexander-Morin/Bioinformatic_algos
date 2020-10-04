# Chapter 1 is centered around text matching algorithms, framed by asking how to find the ori in bacteria.

# Begin with simple string match/count, sliding a window of length pattern down text.
# ----------------------------------------------------------------------------------------------------------------------

def pattern_count(text, pattern):
    count = 0
    k = len(pattern)
    for i in range(0, len(text) - k-1):
        if text[i:i+k] == pattern:
            count += 1
    return count


def pattern_match(text, pattern):
    index_array = []
    k = len(pattern)
    for i in range(0, len(text) - k-1):
        if text[i:i + k] == pattern:
            index_array.append(i)
    return index_array

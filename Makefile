# ============================================================================
#  Makefile - Bermudean Option Pricing Project
#  Compares European, Bermudean, and American option pricing methods
# ============================================================================

CXX      = c++
CXXFLAGS = -std=c++11 -O2 -Wall
TARGET   = projet_bermude

SRCS     = main.cpp \
           BlackScholes.cpp \
           Random.cpp \
           BermudeNDatesLS.cpp \
           PricerArbreBinomial.cpp

OBJS     = $(SRCS:.cpp=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^

%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean

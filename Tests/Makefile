CC      = cc
CFLAGS  = -std=c11 -O2 -I../Sources/PitchDSP/include \
          -DDEBUG_PITCH=0 \
          -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable
SRC     = ../Sources/PitchDSP/PitchDSP.c test_pitchdsp.c
TARGET  = test_pitchdsp

.PHONY: all run clean

all: $(TARGET)

run: $(TARGET)
	./$(TARGET) Resources

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) -lm

clean:
	rm -f $(TARGET)

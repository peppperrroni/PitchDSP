CC      = cc
CFLAGS  = -std=c11 -O2 -I../Sources/PitchDSP/include \
          -DDEBUG_PITCH=0 \
          -Wall -Wno-unused-function -Wno-unused-variable -Wno-unused-but-set-variable
SRC     = ../Sources/PitchDSP/PitchDSP.c test_pitchdsp.c
TARGET  = test_pitchdsp
INTERNAL_TARGET = test_internal

.PHONY: all run clean

all: $(TARGET) $(INTERNAL_TARGET)

run: $(TARGET) $(INTERNAL_TARGET)
	./$(INTERNAL_TARGET)
	./$(TARGET) Resources

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRC) -lm

# test_internal.c #includes PitchDSP.c directly (white-box) — do not link it twice.
$(INTERNAL_TARGET): test_internal.c ../Sources/PitchDSP/PitchDSP.c
	$(CC) $(CFLAGS) -o $(INTERNAL_TARGET) test_internal.c -lm

clean:
	rm -f $(TARGET) $(INTERNAL_TARGET)

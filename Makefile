CC = gcc
CFLAGS = -O3 -Wall -mcpu=native -flto -ffast-math
LDFLAGS = -lm -lpthread -flto

SRCS = spr_fit_main.c spr_io.c spr_models.c spr_optim.c
HDRS = spr_types.h spr_io.h spr_models.h spr_optim.h
TARGET = spr_fit

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SRCS) $(HDRS)
	$(CC) $(CFLAGS) -o $(TARGET) $(SRCS) $(LDFLAGS)

clean:
	rm -f $(TARGET)

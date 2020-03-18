################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
O_SRCS += \
../Calculs.o \
../Coarse.o \
../Compress.o \
../Dump.o \
../Filter.o \
../Main.o \
../Ssh.o \
../Statistics.o \
../Step.o \
../Surface.o 

CPP_SRCS += \
../Calculs.cpp \
../Coarse.cpp \
../Compress.cpp \
../Dump.cpp \
../Filter.cpp \
../Main.cpp \
../Ssh.cpp \
../Statistics.cpp \
../Step.cpp \
../Surface.cpp 

OBJS += \
./Calculs.o \
./Coarse.o \
./Compress.o \
./Dump.o \
./Filter.o \
./Main.o \
./Ssh.o \
./Statistics.o \
./Step.o \
./Surface.o 

CPP_DEPS += \
./Calculs.d \
./Coarse.d \
./Compress.d \
./Dump.d \
./Filter.d \
./Main.d \
./Ssh.d \
./Statistics.d \
./Step.d \
./Surface.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

Coarse.o: ../Coarse.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -O3 -g3 -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"Coarse.d" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '



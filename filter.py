from scipy.io import wavfile
from scipy import signal
import numpy as np
from matplotlib import pyplot as plt
from scipy import interpolate as itp



echo_file_input = "ILoveDSPecho.wav"
tone_file_input = "ILoveDSPtone.wav"

echo_file_output = "FilteredILoveDSPecho.wav"
tone_file_output = "FilteredILoveDSPtone.wav"

def autoCorrelation(timeSignal):
    distanceFromMean = timeSignal - np.mean(timeSignal)
    fftRes = np.fft.fft(distanceFromMean)
    mag = [np.linalg.norm(x)**2 for x in fftRes]
    ifftResReal = np.real(np.fft.ifft(mag))
    return ifftResReal/np.sum(distanceFromMean**2)

def getReal(x):
    return  np.int16(round(x.real))

def outputFormat(x):
    return np.int16(round(x))

def compareSpectrums(data, newData, timeStep, name):
    newData = np.array(list(map(outputFormat, newData)))
    dataFFT = np.fft.fft(data)
    newDataFFT = np.fft.fft(newData)

    dataFrequencies = np.fft.fftfreq(data.size, timeStep)
    newDataFrequencies = np.fft.fftfreq(newData.size, timeStep)
    indices = np.argsort(dataFrequencies)
    newIdicies = np.argsort(newDataFrequencies)

    plt.title("Frequency Analysis for Original Signal and Filtered Signal" + name)
    plt.plot(dataFrequencies[indices], dataFFT[indices], label="Original Signal")
    plt.plot(newDataFrequencies[newIdicies], newDataFFT[newIdicies], label="Filtered Signal")
    plt.xlabel("Frequency")
    plt.ylabel("Signal Strength")
    plt.legend()
    plt.show()

    pass

def remove_echo(file_name, output_file_name):
    print("Reading file for echo removal:\t", file_name)
    fs, data = wavfile.read(file_name)
    

    corr = autoCorrelation(data)

    #second peak sample number
    targetSampleNumber = 22050
    correlationRatio = 0.49165

    #removing echo
    newData = np.array(data)
    for i in range(targetSampleNumber, len(data)):
        newData[i] = data[i] - data[i-targetSampleNumber] * correlationRatio

    newData = np.array(list(map(outputFormat, newData)))


    #recalculating auto correlation
    corrAfterFilter = autoCorrelation(newData)
    
    plt.title("Auto Relation Comparison for original signal and filtered signal")
    plt.plot(corr[:len(corr)//2], label="Original Signal")    
    plt.plot(corrAfterFilter[:len(corr)//2], label="Filtered Signal")
    plt.xlabel("Sample Number")
    plt.ylabel("Correlation Ratio")
    plt.legend()
    plt.show()


    compareSpectrums(data, newData, 1/fs, " after Echo Removal")

    wavfile.write(output_file_name, fs, newData)
    print("Done with writing the file")
    
    pass

def remove_noise(input_file_name, output_file_name):
    
    print("Reading file for noise cancellation:\t", input_file_name)
    fs, data = wavfile.read(input_file_name)

    #same
    newData = data

    #not good frequency
    targetFrequency = 3200
    terminationFrequency = 5000
    
    
    #using butterworth
    #targetMagnitude = 48
    #filter_ = Butter(btype="lowpass", cutoff=int(targetFrequency), rolloff=targetMagnitude, sampling=fs)
    #newData = np.array(filter_.send(list(data)))

    #newData = data - np.mean(data)

    
    #no changes
    #newData = np.fft.ifft(np.fft.fft(data))
    #newData =np.array(list(map(getReal,newData)))


    #using special filter

    #window = np.bartlett(3)
    #newData = np.convolve(newData, window)
    
    timeStep = 1 / fs

    newDataFFT = np.fft.fft(newData)
    newDataFrequencies = np.fft.fftfreq(newData.size, timeStep)
    newIdicies = np.argsort(newDataFrequencies)

    for index in newIdicies:
        if newDataFrequencies[index] > targetFrequency:
            if np.linalg.norm(newDataFFT[index])**2 > (1000000)**2:
                newDataFFT[index] = 0.0000001 * newDataFFT[index]
        elif newDataFrequencies[index] < -targetFrequency:
            if np.linalg.norm(newDataFFT[index])**2 > (1000000)**2:
                newDataFFT[index] = 0.0000001 * newDataFFT[index]

        if newDataFrequencies[index] > terminationFrequency:
            newDataFFT[index] = 0.0 * newDataFFT[index]
        elif newDataFrequencies[index] < -terminationFrequency:
            newDataFFT[index] = 0.0 * newDataFFT[index]

    newData = np.array(list(map(getReal,np.fft.ifft(newDataFFT))))
    

    #using numpy
    #window = np.bartlett(3)
    #newData = np.convolve(newData, window)
    #b, a = signal.butter(1, normalizedFrequency, 'low')
    #newData = signal.lfilter(b, a, newData)
    
    '''
    targetFrequency = 3500
    normalizedFrequency = targetFrequency / (fs//2)
    b, a = signal.butter(10, normalizedFrequency, 'low')
    newData = signal.lfilter(b, a, newData)
    '''

    '''
    targetFrequency = 5500
    normalizedFrequency = targetFrequency / (fs//2)
    b, a = signal.butter(20, normalizedFrequency, 'low')
    newData = signal.lfilter(b, a, newData)
    '''

    #print(len(b))
    #print(b)
    #print(len(a))
    #print(a)

    #smothing newData
    #window = 6 * [1/6]
    #window = np.bartlett(3)
    #window = np.kaiser(3,705)
    #window = np.bartlett(3)
    #newData = np.convolve(newData, window)
    #newData = np.convolve(newData, window)
    #newData = np.convolve(newData, window)
    


    #calcaluting final spectrum for written data
    newData = np.array(list(map(outputFormat, newData)))
    compareSpectrums(data, newData,1/fs, " after noise cancellation")

    
    wavfile.write(output_file_name, fs, newData)
    pass


def filterFile(input_file_name, output_file_name, removeEcho = False, removeNoise = False):
    if removeEcho and removeNoise:
        remove_echo(input_file_name, output_file_name)
        remove_noise(output_file_name, output_file_name)
    elif removeEcho:
        remove_echo(input_file_name, output_file_name)
    elif removeNoise:
        remove_noise(input_file_name, output_file_name)
    pass


def main():
    filterFile(echo_file_input, echo_file_output, True, False)
    filterFile(tone_file_input, tone_file_output, False,  True)


    pass


main()
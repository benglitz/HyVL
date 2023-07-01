# HyVL - Hybrid Vocalization Localizer
The HyVL system relies on specific hardware and software for achieving high accuracy sound localization. We provide guides here for the hardware setup, documentation of the core functions and a sample data analysis pipeline. The present guide assumes that different labs will use different pipelines, programming languages and hardware components, and thus tries to describe general guidelines, rather than minutely reproduce our setup (which is documented in the publication). 

![](SystemOverview.png)

# Hardware Setup
The hardware setup can generally be divided into two components, the high-fidelity microphone setup, the Cam64 setup and the video camera setup. While we provide some details regarding our choice of hardware, the specific devices are of course exchangeable and can be adapted to the specific lab. 
It is generally preferable to keep the whole setup inside a sound attenuated booth with sound absorbing foam on the inside, to avoid echos/reflections of the sound and outside sounds from contaminating the recording.
Since multiple data sources are being acquired with substantial data throughput each, which puts the computer under considerable load, it is recommendable to not use the computer for other tasks during the acquisition process.
Generally, all measurements are supplied in units of meters to avoid conversion problems.

## High-Fidelity Microphone Setup
An array of high-fidelity ultrasonic microphones forms the basis for the acquisition of the audio data for USV detection and SLIM sound localization. While the geometry of the array and the number of microphones can be chosen rather freely (our code in VocLocalizer.m adapts to these), the accuracy of the design is essential for high-precision localization. 
- Number of Microphones: 3 is the required minimum, but the more the better. We used 4 microphones, but expect that each additional microphone will increase the spatial resolution further.
- Type of Microphone: We used Avisoft CMPA/16 microphones due to their high degree of spectral range and flatness, high sensitivity and possibility to directly output to DAQ hardware. Using USB-based microphones is also possible but would require additional synchronization, similar to the synchronization for Cam64 (see below). 
- Sampling Rate: 250 kHz or more should be sufficient in mice and rodents.
- Positioning: The positioning of the microphones can be customized by the user and tailored to the recording environment. The specific positions are then passed to the VocLocalizer, which takes the specific positions into account.

For positioning, the main limitations to take into account are: 
  - The locations of the microphones should encompass the recording area to maximize runtime differences.
   -The microphones should be placed higher than the platform to avoid sound absorption by the bodies of the mice.
  - Ultrasound is rather directional, therefore the microphones should be distributed in a way to maximize the chance of being in the direction of vocalization of the mice.

## Cam64 Setup
The Cam64 performs the acquisition of the high-channel count audio data which is then used for beamforming, on the basis of the detected USV times and frequency ranges.
- Sampling Rate: 250 kHz is the maximal SR supported by the Cam64 and an excellent fit for the recordings of mouse USVs.
- Positioning: The Cam64 should be installed directly above the platform. For simplicity the base of the Cam64 should be parallel to the surface of the interaction platform. Importantly, minimize all deviations/rotations of the Cam64, which will otherwise impact accuracy later. The distance of the Cam64 to the platform should be as small as possible to improve the SNR or the vocalizations. We used a distance of 46.5cm, but a small distance would likely have resulted in better SNR of the USVs.
- Cam64 internal camera: The camera of the Cam64 is in principle able to provide images for later video tracking. We, however, noticed that the camera acquisition creates a high-frequency noise, in the range of about 50kHz, and thus interferes with audio recording of USVs. In addition, the frames from the Cam64 camera are not precisely timed, and thus will only be approximately aligned to the audio data. We therefore chose to use an external camera instead.

## Video Camera Setup
The Video camera setup can be customized to the specific models/systems available in the lab. We used a PointGrey monochrome camera, which was supported by the image acquisition toolbox in Matlab and allowed for triggering of individual frames.
- Positioning: The camera should be placed with the line of sight orthogonal to the platform, in order to avoid any spatial non-uniformities/scalings in different locations of the image.
- Lighting: In our experiments we used white lighting which seemed to not irritate the mice, however, some researchers might prefer (infra)red light, which is less or invisible to the mice and thus interferes less with their natural nocturnal behavior. White light can provide advantages for tracking using a color camera (see below).
- Camera settings: To simplify later tracking, we would recommend to use high resolution acquisition from the start. In hindsight, our resolution (640x512) was likely on the low end and tracking would improve with higher resolution. Also, using a color camera will provide additional information in the different channels, which can be beneficial for tracking in particular during close interaction. We used a framerate of 50 fps, which appears to have been a good choice. Exposure duration should be kept short in order to avoid blurring the images.

## Integration and Synchronization
Synchronization between the different data sources is essential, given that the behavior or the mice is rather rapid and the duration of the USVs rather short, down to 10 ms. Since the different data sources required different synchronization strategies, we detail them separately here.

### Temporal Synchronization
Overall we recommend using a single data acquisition device that maintains a global, precise time, which all other devices are then synchronized to. We used an NI DAQ card (NI-6321), which acquired synchronization pulses at 10 kHz and provided trigger pulses to other devices . 
- USM4 Data: We used a second NI DAQ (NI-6259) card to acquire 4 channels of audio at 250 kHz. Both cards were synched internally and linked through a trigger cable. The acquisition process should be triggered/aligned precisely to the synchronization process on the first card.
- Cam64 Data: The Cam64 system does not support precisely times/triggered starting. As an alternative one therefore needs to provide an acoustic synchronization pulse. This can be achieved by placing a miniature speaker (e.g. old in-ear headphone driver) next to one of the microphones and providing a very brief sound (e.g. rapid voltage step). The speaker can be directly connected to an digital output on the DAQ card, without the need of an additional amplifier. We provide a spatial mapping of the Cam64 microphones here, for choosing and later identifying a suitable microphone. We chose microphone 58 on the outside of the array.
- Video Data: The video camera acquired a single image per TTL pulse it received and also provided a TTL pulse to confirm the acquisition of the device. These pulses were all provided and acquired by the main DAQ-card. This allows post-hoc precise timing of each frame. Importantly, frame loss can sometimes occur, which is not visible on the level of these triggers, e.g. if the USB interface drops frames. The driver of the Image Acquisition Toolbox provides the timing of the received frames, which can be compared to the triggered frames to identify the missing frames.

### Spatial Synchronization
We converted all estimates to a general coordinate system whose origin was placed in the center of the interaction platform. This simplified later integration of all spatial estimates.
- USM4 Data: The estimates between every pair of microphones generates an estimate on line connecting the microphones, specifically on the height of the platform. These estimates are combined in VocCollector and based on the microphone positions (which are in the coordinates mentioned above) to a single estimate on the platform (or rather, elevated by about 1cm to the height of the snout of the mice).
- Cam64 Data: The estimates of the Cam64 are simply shifted by the position of the Cam64 relative to the platform to obtain them in the coordinates mentioned above.
- Video Data: An important step in this transition from pixel-coordinates to flat, Euclidean coordinates is to unwarp any lens-distortions of the acquired image. In our case, the camera/lens was already relatively flat and only mild unwarping was required. We provide multiple functions for this purpose in the repository (`C_correctImageDistortion`, `C_correctImageDistortion_V1_1`).




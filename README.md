# Vision Based Dynamic Offside Line Marker for Soccer Games

Link to detailed description: [Vision Based Dynamic Offside Line Marker for Soccer Games](https://arxiv.org/pdf/1804.06438.pdf)

Offside detection in soccer has emerged as one
of the most important decision with an average of 50 offside
decisions every game. False detections and rash calls adversely
affect game conditions and in many cases drastically change
the outcome of the game. The human eye has finite precision
and can only discern a limited amount of detail in a given
instance. Current offside decisions are made manually by
sideline referees and tend to remain controversial in many
games. This calls for automated offside detection techniques
in order to assist accurate refereeing.


In this work, we have explicitly used computer vision and
image processing techniques like Hough transform, color similarity (quantization), graph
connected components, and vanishing point ideas to identify
the probable offside regions.

Due to the unavailability of specific datasets, we have implemented and evaluated our algorithms on FIFA 16 game videos with a standard camera setting.

### Algorithmic Workflow
![Dropouts](https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/Workflow.png)

### Identifying the Field of Play using Hough Transform and Colour Similarity
| Input Video   | Cropped Field of Play|
| ------------- |:-------------:|
|  <img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/Hough_Input.gif" width="400">   | <img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/Hough_Output.gif" width="400"> |
|  <img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/Input_Aud.PNG" width="400">   | <img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/Output_Aud.PNG" width="400"> |
### Localizing Players within the Field of Play
| Input Video | Identifying Attacking Team  | Localizing player locations|
| ------------- |:-------------:|:-------------:|
|  <img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/image10.gif" width="400">   | <img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/image14.gif" width="400">   | <img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/image15.gif" width="400"> |
### Tracking the Last Defender using the KLT tracker
###### Figure shows the efficient team-wise player classification & tracking of the last defender to determine the offside marker

<img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/image20.gif" width="800">  

### Marking the Offside Line - Vanishing Points

| Input Video   | Cropped Field of Play|
| ------------- |:-------------:|
|  <img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/image17.png" width="400">   | <img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/image18.png" width="400"> |

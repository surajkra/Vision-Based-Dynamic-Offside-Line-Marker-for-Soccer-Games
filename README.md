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



| Input Video   | Cropped Field of Play|
| ------------- |:-------------:|
|  !<img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/Hough_Input.gif" width="400">   | <img src="https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/Hough_Output.gif" width="400"> |

| Tables        | Are           | Cool  |
| ------------- |:-------------:| -----:|
| col 3 is      | right-aligned | $1600 |
| col 2 is      | centered      |   $12 |
| zebra stripes | are neat      |    $1 |

### Identifying the Field Boundaries using Hough Transform
![Dropouts](https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/Hough_Input.gif){:height="50%" width="50%"}
![Dropouts](https://github.com/surajkra/Vision-Based-Dynamic-Offside-Line-Marker-for-Soccer-Games/blob/master/Images/Hough_Output.gif)

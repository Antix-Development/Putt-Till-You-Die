<div align="justify">

# Putt Till You Die

<p align="center"><img src="https://github.com/Antix-Development/Putt-Till-You-Die/blob/master/img/450X250.png"></p>

<br>

Putt Till You Die is my first submission for the 2022 [js13kGames](https://js13kgames.com) competition, where over a one month period you create a game based on a theme (dEATH FOR 2022) and then you cram your entire game (code and assets) into a 13Kb zip file and submit it.

You can play Putt Till You Die on [its js13k competition page](https://js13kgames.com/entries/), or download and load index.html from this repository.

<br>

## Introduction

Putt Till You Die is a 2 dimensional bird's eye view mini-putt golf game written in JavaScript, HTML, and CSS.

The game its self doesn't really fit in with the theme but hey, you could literally die whilst playing the game.. as you could any other game ;)

Also during the production of this game, my dear friend *Ewook* passed away after a short battle with cancer. He loved golf and up until a month before his passing was playing 4 rounds a week. He was 84 and this game is dedicated to him.

<br>

## Gameplay

The aim of Putt Till You Die is to basically putt the ball into the goal. That's prety much it.

<br>

## Features

- Unlimited procedural generated holes!
- Seeded PRNG (Pseudo Random Number Generator) so the game is the same across devices.
- A persistent game state.

<br>

## Technologies

Putt TIll You Die was developed using Visual Studio Code and Notepad++.

Particle imagery was created using GIMP, and all other imagery is drawn programatically.

Sound effects were created using Frank Force's fantastic [ZzFX](https://killedbyapixel.github.io/ZzFX/).

[Base64 Image](https://www.base64-image.de) was used to convert the PNG particle image to an embeddable data URI.

[Terser](https://terser.org/) was used to minify the JavaScript code.

Putt Till You Die has been tested working in FireFox and Chrome browsers, your milage might vary with other browsers.

<br>

## The Code

[Mini Physics](https://github.com/xem/mini2Dphysics) by Xem us used for the physics simulation. It has a small issue where if a circle collides with two axis aligned boxes that are directly adjacent to eachother, the collision does not work correctly. If I had more time I'd just write my own code but time does not allow :(

I'm always more than happy to answer any questions you might have regarding the code, just drop me a line at antix.development@gmail.com :)

</div>
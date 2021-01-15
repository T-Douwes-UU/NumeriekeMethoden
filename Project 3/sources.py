"""
This file stores functions that are used across multiple tasks, in order of usage.
"""
from matplotlib.animation import FuncAnimation
import mpl_toolkits.axes_grid1
import matplotlib.widgets


class Player(FuncAnimation):
    """Matplotlib video player, adapted from code courtesy of Elan Ernest."""
    def __init__(self, fig, func, frames=None, start=0, init_func=None, fargs=None,
                 save_count=None, pos=(0.125, 0.92), **kwargs):
        self.min = 0
        self.scale = 1
        if isinstance(frames, int):
            self.max = frames - 1
        elif isinstance(frames, range):
            self.min = frames.start
            self.max = frames.stop - 1
            self.scale = frames.step
        else:
            self.max = 100
        self.i = start if self.min <= start <= self.max else self.min
        self.step = self.scale
        self.fig = fig
        self.func = func
        self.artists = self.func(self.i)
        self.setup(pos)
        FuncAnimation.__init__(self, self.fig, self.update, frames=self.play(),
                               init_func=init_func, fargs=fargs,
                               save_count=save_count, **kwargs)

    def play(self):
        while self.step:
            self.i += self.step
            if self.i < self.min:
                self.i = self.max
            elif self.i > self.max:
                self.i = self.min
            yield self.i

    def start(self):
        if not self.step:
            self.event_source.start()

    def stop(self, event=None):
        if self.step:
            self.step = 0
            self.event_source.stop()

    def forward(self, event=None):
        if self.step > self.scale:
            self.step //= 2
        else:
            self.start()
            self.step = self.scale

    def backward(self, event=None):
        if self.step < -self.scale:
            self.step //= 2
        else:
            self.start()
            self.step = -self.scale

    def fastforward(self, event=None):
        if self.step > 0:
            self.step *= 2
        else:
            self.start()
            self.step = 2 * self.scale

    def fastback(self, event=None):
        if self.step < 0:
            self.step *= 2
        else:
            self.start()
            self.step = -2 * self.scale

    def oneforward(self, event=None):
        if self.i < self.max:
            self.i += self.scale
        self.onestep()

    def onebackward(self, event=None):
        if self.i > self.min:
            self.i -= self.scale
        self.onestep()

    def onestep(self):
        self.stop()
        self.func(self.i)
        self.slider.set_val(self.i)
        self.fig.canvas.draw_idle()

    def setup(self, pos):
        playerax = self.fig.add_axes([pos[0], pos[1], 0.64, 0.04])
        divider = mpl_toolkits.axes_grid1.make_axes_locatable(playerax)
        fbax = divider.append_axes("right", size="80%", pad=0.03)
        bax = divider.append_axes("right", size="80%", pad=0.03)
        sax = divider.append_axes("right", size="80%", pad=0.03)
        fax = divider.append_axes("right", size="80%", pad=0.03)
        ffax = divider.append_axes("right", size="80%", pad=0.03)
        ofax = divider.append_axes("right", size="100%", pad=0.03)
        sliderax = divider.append_axes("right", size="400%", pad=0.05)
        self.button_oneback = matplotlib.widgets.Button(playerax, label='$\u29CF$')
        self.button_fastback = matplotlib.widgets.Button(fbax, label='$\u25C2\u25C2$')
        self.button_back = matplotlib.widgets.Button(bax, label='$\u25C0$')
        self.button_stop = matplotlib.widgets.Button(sax, label='$\u25A0$')
        self.button_forward = matplotlib.widgets.Button(fax, label='$\u25B6$')
        self.button_fastforward = matplotlib.widgets.Button(ffax, label='$\u25B8\u25B8$')
        self.button_oneforward = matplotlib.widgets.Button(ofax, label='$\u29D0$')
        self.button_oneback.on_clicked(self.onebackward)
        self.button_fastback.on_clicked(self.fastback)
        self.button_back.on_clicked(self.backward)
        self.button_stop.on_clicked(self.stop)
        self.button_forward.on_clicked(self.forward)
        self.button_fastforward.on_clicked(self.fastforward)
        self.button_oneforward.on_clicked(self.oneforward)
        self.slider = matplotlib.widgets.Slider(sliderax, '',
                                                self.min, self.max, valinit=self.i)
        self.slider.on_changed(self.set_pos)

    def set_pos(self, i):
        self.i = int(self.slider.val)
        self.artists = self.func(self.i)

    def update(self, i):
        self.slider.set_val(i)
        return self.artists
